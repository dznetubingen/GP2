#!/usr/bin/env python
import glob
import os
import pandas as pd
from pathlib import Path
from snakemake.utils import min_version
from collections import defaultdict

shell.prefix("set -o pipefail; umask 002; ")  # set g+w

# create directories that we need to store the outputs
DIRs = ['vcfs_per_fam', 'short_list', 'reports', 'utr']

for DIR in DIRs:
    if not os.path.exists(DIR):
        os.makedirs(DIR)

configfile: "configs/config.yaml"

###make tfam file for each family
df= pd.read_csv(config["gp2_ped"], sep='\t', header=None, names=['family_id','id','parental_id','maternal_id','sex','phenotype'])
#wgs = pd.read_csv('inputs/wgs_samples.list', header=None, names=['id'])

cohort=defaultdict(list,{ k:[] for k in ('singleton','duo','group','trio','multi') })
for family, dat in df.groupby('family_id'):
    if (dat['phenotype'] == 1).all():
        print(f"No affected samples in {family}.")
    else:
        if dat.shape[0] == 1:
            cohort['singleton'].append(family)
        elif dat.shape[0] == 2:
        ## this is to deal with the other unaffected family member who is not a parent of the affected sample
            if (dat.iloc[:,2:4]== "0").all(axis=None) and dat.iloc[:,-1].sum() == 3:
                cohort['group'].append(family)
            elif dat.iloc[:,-1].sum() == 4:
                 cohort['multi'].append(family)
            else:
                cohort['duo'].append(family)
        elif dat.shape[0] == 3 and dat['phenotype'].sum() ==4:
            cohort['trio'].append(family)
        else:
            cohort['multi'].append(family) 

#for fam in list(cohort.keys()): 
#    N=len(cohort[fam])
#    print(f"N of {fam} in cohort:{N}")

ped_dict = {v: k for k, values in cohort.items() for v in values}

localrules: all, merge_tsv, make_tfam

ruleorder: make_tfam > subset_per_fam


rule make_tfam:
    output:
        tfam= ["vcfs_per_fam/{famID}.tfam".format(famID=famID) for famID in ped_dict.keys()],
        fam_list= ["vcfs_per_fam/{famID}.list".format(famID=famID) for famID in ped_dict.keys()] 
    run:
        import pandas as pd
        df= pd.read_csv(config["gp2_ped"], sep="\t", header=None, names=["family_id","id","parental_id","maternal_id","sex","phenotype"])
        FAM = list(df['family_id'].unique())

        for family, dataframe in df.groupby('family_id'):
            # save the dataframe for each group to a csv
            dataframe.to_csv(f'vcfs_per_fam/{family}.tfam', sep="\t", header=False, index=False)
            dataframe[['id']].to_csv(f'vcfs_per_fam/{family}.list', sep="\t", header=False, index=False)

rule subset_per_fam:
    input: 
        vcf = "final_vcfs/GP2_annotate.vcf.gz",
        fam_list = "vcfs_per_fam/{famID}.list"
    output:
        "vcfs_per_fam/raw/{famID}.vcf.gz"
    conda:
        "../envs/tools.yml"
    threads: 1
    resources:  
        nodes = 1
    shell:
        """
        ##subset the vcf by family and remove monomorphic ref sites in the family vcf
        bcftools view -S {input.fam_list} {input.vcf} \
        |  bcftools filter -e 'AC==0' -O z -o {output} && tabix -p vcf -f {output}
        """

rule slivar_filter:
    input:
        vcf = rules.subset_per_fam.output,
        ped = "vcfs_per_fam/{famID}.tfam"
    output:
        temp("vcfs_per_fam/seg/{famID}.vcf")
    threads: 1
    resources:
        nodes = 1
    run:
        ## trio: call recessive variants for the complete penetrance and call dominant for incomplete penetrance 
        if ped_dict[wildcards.famID] == 'trio':
            shell("""singularity exec -B /mnt/vol009/GP2/vep:$HOME -B /mnt/vol009/hg38_vep105:/mnt slivar_modified_0.2.7.sif \
                  slivar expr --vcf {input.vcf} \
                    --ped {input.ped} \
                    -o {output} \
                    --js /mnt/slivar/slivar-functions.js \
                    --pass-only \
                    --info "variant.ALT[0] != '*' && variant.call_rate >= 0.95" \
                    --trio "denovo:(variant.CHROM != 'X' || variant.CHROM != 'chrX') && denovo(kid, mom, dad) && INFO.gnomad_popmax_af < 0.01" \
                    --trio "x_denovo:x_denovo(kid, mom, dad) && INFO.gnomad_popmax_af <= 0.01" \
                    --trio "HOM:(variant.CHROM != 'X' || variant.CHROM != 'chrX') && recessive(kid, mom, dad) && INFO.gnomad_nhomalt < 10" \
                    --trio "X_HOM:x_recessive(kid, mom, dad) && INFO.gnomad_nhomalt <= 10" \
                    --trio "comphet:comphet_side(kid, mom, dad) && INFO.gnomad_nhomalt <= 10" \
                    --family-expr "dominant:(variant.CHROM != 'X' || variant.CHROM != 'chrX') && fam.every(segregating_dominant) && !fam.every(segregating_denovo) && INFO.gnomad_popmax_af <= 0.01" \
                    --family-expr "x_dominant:(variant.CHROM == 'X' || variant.CHROM == 'chrX') && fam.every(segregating_dominant_x) && !fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af <= 0.01"
                  """)
        else: 
            shell("""singularity exec -B /mnt/vol009/GP2/vep:$HOME -B /mnt/vol009/hg38_vep105:/mnt slivar_modified_0.2.7.sif \
                  slivar expr --vcf {input.vcf} \
                    --ped {input.ped} \
                    -o {output} \
                    --js /mnt/slivar/slivar-functions.js \
                    --pass-only \
                    --info "variant.ALT[0] != '*' && variant.call_rate >= 0.95" \
                    --family-expr "HOM:(variant.CHROM != 'X' || variant.CHROM != 'chrX') && fam.every(segregating_recessive) && INFO.gnomad_nhomalt <= 10" \
                    --family-expr "X_HOM:(variant.CHROM == 'X' || variant.CHROM == 'chrX') && fam.every(segregating_recessive_x) && INFO.gnomad_nhomalt <= 10" \
                    --family-expr "denovo:(variant.CHROM != 'X' || variant.CHROM != 'chrX') && fam.every(segregating_denovo) && INFO.gnomad_popmax_af <= 0.01" \
                    --family-expr "x_denovo:(variant.CHROM == 'X' || variant.CHROM == 'chrX') && fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af <= 0.01" \
                    --family-expr "dominant:(variant.CHROM != 'X' || variant.CHROM != 'chrX') && fam.every(segregating_dominant) && !fam.every(segregating_denovo) && INFO.gnomad_popmax_af <= 0.01" \
                    --family-expr "x_dominant:(variant.CHROM == 'X' || variant.CHROM == 'chrX') && fam.every(segregating_dominant_x) && !fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af <= 0.01" \
                    --family-expr "comphet_side:fam.every(function(s) {{return (s.het || s.hom_ref) && hq1(s) }}) && fam.some(function(s) {{return s.het && s.affected}}) && INFO.gnomad_nhomalt <= 10"
                   """)

rule Snpsift_count:
    input:
        vcf = rules.slivar_filter.output,
        tfam = "vcfs_per_fam/{famID}.tfam"
    output:
        vcf = temp("vcfs_per_fam/count/seg_{famID}.vcf")
    conda:
        "../envs/tools.yml"
    threads: 1
    resources:
        nodes = 1
    shell:
        """
        SnpSift caseControl -tfam {input.tfam} {input.vcf} > {output.vcf}
        """

skip_list = [
    'intergenic',
    '5_prime_UTR',
    '3_prime_UTR',
    'TFBS_ablation',
    'TF_binding_site',
    'regulatory_region',
    'non_coding_transcript',
    'non_coding',
    'upstream_gene',
    'downstream_gene',
    'non_coding_transcript_exon',
    'NMD_transcript'
    ]

rule slivar_comphet:
    input:
        vcf = rules.Snpsift_count.output.vcf,
        ped = "vcfs_per_fam/{famID}.tfam"
    output:
        "vcfs_per_fam/comphet/{famID}.vcf"
    params:
        skip = ",".join(skip_list)
    container:
        "docker://zihhuafang/slivar_modified:0.2.7"
    threads: 1
    resources:
        nodes = 1
    shell:
        """
        slivar compound-hets \
               --skip {params.skip} \
               -v {input.vcf} \
               --allow-non-trios --sample-field comphet_side --sample-field denovo -p {input.ped} \
               -o {output}
        """

info_fields = [
    'gnomad_AF',
    'gnomad_popmax_af',
    'gnomad_nhomalt',
    'gnomad_AC',
    'TOPMed8_AF',
    'CADD_PHRED',
    'GP2_affected_gt',
    'amp_pd_cases_gt',
    'Cases',
    'Controls',
    'CNCR'
    ]

csq_columns= [
    'Existing_variation',
    'MAX_AF_POPS',
    'LoF',
    'LoF_filter',
    'LoF_flags',
    'LoF_info',
    'MetaRNN_score',
    'Ensembl_transcriptid',
    'SpliceAI_pred_DS_AG',
    'SpliceAI_pred_DS_AL',
    'SpliceAI_pred_DS_DG',
    'SpliceAI_pred_DS_DL',
    'SpliceAI_pred_SYMBOL',
    'SpliceRegion',
    'existing_InFrame_oORFs',
    'existing_OutOfFrame_oORFs',
    'existing_uORFs',
    'five_prime_UTR_variant_annotation',
    'five_prime_UTR_variant_consequence',
    'Gene',
    'IMPACT',
    'BIOTYPE',
    'STRAND',
    'CANONICAL',
    'MANE_SELECT',
    'MANE_PLUS_CLINICAL',
    'TSL',
    'NEAREST',
    'EXON',
    'Codons',
    'Amino_acids',
    'HGVSc',
    'HGVSp',
    'ClinVar_CLNSIG',
    'ClinVar_CLNDN',
    'VAR_SYNONYMS',
    'MOTIF_NAME',
    'MOTIF_POS',
    'HIGH_INF_POS',
    'MOTIF_SCORE_CHANGE',
    'TSSDistance',
    'miRNA',
    'TRANSCRIPTION_FACTORS',
    'NMD',
    'DisGeNET',
    'DownstreamProtein',
    'ProteinLengthChange'
    ]

gene_descs= [
    '/mnt/slivar/pLI_lookup.txt',
    '/mnt/slivar/oe_lof_upper_lookup.txt',
    '/mnt/slivar/oe_mis_upper_lookup.txt',
    '/mnt/slivar/oe_syn_upper_lookup.txt',
    '/mnt/slivar/clinvar_gene_desc.txt',
    '/mnt/slivar/hgncSymbol.inheritance.tsv',
    '/mnt/slivar/lookup_description_symbol.txt',
    '/mnt/slivar/OMIM_entry_genesymbol.txt'
    ]

rule slivar_tsv:
    input:
        vcf = rules.Snpsift_count.output.vcf,
        ped = "vcfs_per_fam/{famID}.tfam"
    output:
        "short_list/filtered_{famID}.tsv"
    params:
        info = "".join([f"--info-field {x} " for x in info_fields]),
        csq = "".join([f"--csq-column {x} " for x in csq_columns]),
        gene_desc = "".join([f"-g {x} " for x in gene_descs])
    container:
        "docker://zihhuafang/slivar_modified:0.2.7"
    threads: 1
    resources:
        nodes = 1
    shell:
        """
        slivar tsv \
               {params.info} \
               -s denovo \
               -s x_denovo \
               -s HOM \
               -s X_HOM \
               -s dominant \
               -s x_dominant \
               -c CSQ \
               {params.csq} \
               {params.gene_desc} \
               -p {input.ped} \
               {input.vcf} > {output}        
        """

rule comphet_tsv:
    input:
        vcf = rules.slivar_comphet.output,
        ped = "vcfs_per_fam/{famID}.tfam"
    output:
        "short_list/comphet_{famID}.tsv"
    params:
        info = "".join([f"--info-field {x} " for x in info_fields]),
        csq = "".join([f"--csq-column {x} " for x in csq_columns]),
        gene_desc = "".join([f"-g {x} " for x in gene_descs])
    container:
        "docker://zihhuafang/slivar_modified:0.2.7"
    threads: 1
    resources:
        nodes = 1
    shell:
        """
        slivar tsv \
            -s slivar_comphet \
               {params.info} \
               -c CSQ \
               {params.csq} \
               {params.gene_desc} \
             -p {input.ped} \
             {input.vcf} \
             | {{ grep -v ^# || true; }} >> {output}
        """

rule merge_tsv:
    input: 
        filtered = "short_list/filtered_{famID}.tsv",
        comphet = "short_list/comphet_{famID}.tsv"
    output:
        "short_list/all_{famID}.tsv"
    threads: 1
    resources:
        nodes = 1
    shell:
        """
        # get header from first file and drop it from other files
        # and make sure slivar_comphet id is unique
        #awk 'NR == FNR || FNR > 1 {{ sub(/^slivar_comphet/, "comphet", $0); print; }}' {input.filtered} {input.comphet} 
        ## change slivar_comphet_* to comphet 
        cat {input.filtered} {input.comphet} | sed -r -e 's/slivar_comphet_[0-9]+/comphet/g' \
        | sed '1 s/gene_description_1/gnomAD_pLI/;s/gene_description_2/gnomAD_oe_lof_CI90/;s/gene_description_3/gnomAD_oe_mis_CI90/;s/gene_description_4/gnomAD_oe_syn_CI90/;s/gene_description_5/clinvar_gene_description/;s/gene_description_6/MOI/;s/gene_description_7/gene_fullname/;s/gene_description_8/OMIM_link/;' > {output}
        """

rule report:
    input: 
        tsv = "short_list/all_{famID}.tsv",
        ped = "inputs/gp2_all.ped"
    output:
        excel="reports/{famID}.xlsx"
    threads: 1
    script:
        ### modify the report script
        "../scripts/report_format_v4.py"

rule literature:
    input:
        rules.report.output
    output:
        excel="reports/{famID}_pubmed_table.xlsx"
    conda:
        "../envs/literature.yml"
    threads: 1
    shell:
        """
        Rscript --vanilla scripts/pubmed_search.R {input} {wildcards.famID}
        """
