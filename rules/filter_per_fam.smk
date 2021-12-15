#!/usr/bin/env python
import glob
import os
import pandas as pd
from pathlib import Path
from snakemake.utils import min_version
from collections import defaultdict

shell.prefix("set -o pipefail; umask 002; ")  # set g+w

if not os.path.exists("logs"):
    os.makedirs("logs")

if not os.path.exists("short_list"):
    os.makedirs("short_list")

if not os.path.exists("reports"):
    os.makedirs("reports")


#if not os.path.exists("vcfs_per_fam"):
#     os.makedirs("vcfs_per_fam")

configfile: "configs/config.yaml"

###make tfam file for each family
df= pd.read_csv(config["gp2_ped"], sep="\t", header=None, names=["family_id","id","parental_id","maternal_id","sex","phenotype"])
FAM= list(df.family_id.unique())

cohort=defaultdict(list,{ k:[] for k in ('singleton','duo','trio','multi') })
for family, dataframe in df.groupby('family_id'):
    if dataframe.shape[0] == 1:
       cohort['singleton'].append(family)
    elif dataframe.shape[0] == 2:
        ## if parent ids are all 0 then not duo
        if (dataframe.iloc[:,2:4]== "0").all(axis=None) or dataframe.iloc[:,-1].sum() == 4:
            cohort['multi'].append(family)
        else:
            cohort['duo'].append(family)
    elif dataframe.shape[0] == 3:
        ## if parent ids are all 0 then not trio or parent(s) are affected
        if not (dataframe.iloc[:,2:4]== "0").all(axis=None) and dataframe.iloc[:,-1].sum() ==4:
            cohort['trio'].append(family)
        else:
            cohort['multi'].append(family)
    else:
        cohort['multi'].append(family)

ped_dict = {v: k for k, values in cohort.items() for v in values}

localrules: all, merge_tsv


rule subset_per_fam:
    input: 
        vcf = "final_vcfs/GP2_annotate.vcf.gz",
        fam_list = "vcfs_per_fam/{famID}.list"
    output:
        temp("vcfs_per_fam/raw/{famID}.vcf")
    conda:
        "../envs/tools.yml"
    resources:
        nodes = 1
    shell:
        """
        bcftools view -S {input.fam_list} {input.vcf} | SnpSift filter "(countVariant() > 0)" > {output}
    
        """

rule slivar_filter:
    input:
        vcf = rules.subset_per_fam.output,
        ped = "vcfs_per_fam/{famID}.tfam"
    output:
        temp("vcfs_per_fam/seg/{famID}.vcf")
    resources:
        nodes = 1
    run:
        if ped_dict[wildcards.famID] == 'trio':
            shell("""singularity exec -B /mnt/vol009/GP2/vep:$HOME slivar_v0.2.7.sif \
                  slivar expr --vcf {input.vcf} \
                    --ped {input.ped} \
                    -o {output} \
                    --js slivar/slivar-functions.js \
                    --pass-only \
                    --info 'variant.ALT[0] != "*" && variant.call_rate > 0.95 && INFO.gnomad_af < 0.05' \
                    --trio 'denovo:denovo(kid, mom, dad) && INFO.gnomad_nhomalt < 10' \
                    --trio 'x_denovo:x_denovo(kid, mom, dad) && INFO.gnomad_nhomalt < 10' \
                    --trio 'recessive:recessive(kid, mom, dad) && INFO.gnomad_nhomalt < 10' \
                    --trio 'x_recessive:x_recessive(kid, mom, dad)' && INFO.gnomad_nhomalt < 10 \
                    --trio 'comphet:comphet_side(kid, mom, dad)' && INFO.gnomad_nhomalt < 10
                  """)
        else:
            shell("""singularity exec -B /mnt/vol009/GP2/vep:$HOME slivar_v0.2.7.sif \
                  slivar expr --vcf {input.vcf} \
                    --ped {input.ped} \
                    -o {output} \
                    --js slivar/slivar-functions.js \
                    --pass-only \
                    --info 'variant.ALT[0] != "*" && variant.call_rate > 0.95 && INFO.gnomad_af < 0.05' \
                    --family-expr 'recessive:fam.every(segregating_recessive) && INFO.gnomad_nhomalt < 10' \
                    --family-expr 'x_recessive:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_recessive_x) && INFO.gnomad_nhomalt < 10' \
                    --family-expr 'denovo:fam.every(segregating_denovo) && INFO.gnomad_nhomalt < 10' \
                    --family-expr 'x_denovo:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_denovo_x) && INFO.gnomad_nhomalt < 10' \
                    --family-expr 'dominant:INFO.gnomad_af < 0.01 && fam.every(segregating_dominant) && !fam.every(segregating_denovo)' \
                    --family-expr 'x_dominant:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_dominant_x) && !fam.every(segregating_denovo_x) && INFO.gnomad_af < 0.01' \
                    --family-expr 'comphet_side:fam.every(function(s) {{return (s.het || s.hom_ref) && hq1(s) }}) && fam.some(function(s) {{return s.het && s.affected}}) && INFO.gnomad_nhomalt < 10'
                   """)

rule Snpsift_count:
    input:
        vcf = rules.slivar_filter.output,
        tfam = "vcfs_per_fam/{famID}.tfam"
    output:
        temp("vcfs_per_fam/count/seg_{famID}.vcf")
    conda:
        "../envs/tools.yml"
    resources:
        nodes = 1
    shell:
        """
        SnpSift caseControl -tfam {input.tfam} {input.vcf} > {output}

        """


skip_list = [
    'intron',
    'non_coding_transcript',
    'non_coding',
    'upstream_gene',
    'downstream_gene',
    'non_coding_transcript_exon',
    'NMD_transcript'
    ]


rule slivar_comphet:
    input:
        vcf = rules.Snpsift_count.output,
        ped = "vcfs_per_fam/{famID}.tfam"
    output:
        "vcfs_per_fam/comphet/{famID}.vcf"
    params:
        skip = ",".join(skip_list)
    container:
        "docker://brentp/slivar:v0.2.7"
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
    'gnomad_af',
    'gnomad_nhomalt',
    'gnomad_ac',
    'topmed_af',
    'CADD_RAW',
    'CADD_PHRED',
    'GP2_affected_gt',
    'GP2_unaffected_gt',
    'amp-pd_case_nhet',
    'amp-pd_case_nhomalt',
    'amp-pd_case_nmiss',
    'Cases',
    'Controls'
    ]

csq_columns= [
    'Existing_variation',
    'MAX_AF_POPS',
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
    'LoFtool',
    'Gene',
    'BIOTYPE',
    'STRAND',
    'CANONICAL',
    'NEAREST',
    'EXON',
    'Codons',
    'Amino_acids',
    'HGVSc',
    'HGVSp',
    'clinvar_clnsig',
    ]


rule slivar_tsv:
    input:
        vcf = rules.Snpsift_count.output,
        ped = "vcfs_per_fam/{famID}.tfam"
    output:
        "short_list/filtered_{famID}.tsv"
    params:
        info = "".join([f"--info-field {x} " for x in info_fields]),
        csq = "".join([f"--csq-column {x} " for x in csq_columns])
    container:
        "docker://brentp/slivar:v0.2.7"
    resources:
        nodes = 1
    shell:
        """
        n_sample=$(cat {input.ped} | wc -l)
        if [ $n_sample == "1" ];
        then
            slivar tsv \
               {params.info} \
               -i genic \
               -s recessive \
               -s x_recessive \
               -s dominant \
               -s x_dominant \
               -c CSQ \
               {params.csq} \
               -g slivar/lof_lookup.txt \
               -g slivar/clinvar_gene_desc.txt \
               -p {input.ped} \
               {input.vcf} > {output}
        else
            slivar tsv \
               {params.info} \
               -i genic \
               -s denovo \
               -s x_denovo \
               -s recessive \
               -s x_recessive \
               -s dominant \
               -s x_dominant \
               -c CSQ \
               {params.csq} \
               -g slivar/lof_lookup.txt \
               -g slivar/clinvar_gene_desc.txt \
               -p {input.ped} \
               {input.vcf} > {output}        
        fi
        """

rule comphet_tsv:
    input:
        vcf = rules.slivar_comphet.output,
        ped = "vcfs_per_fam/{famID}.tfam"
    output:
        "short_list/comphet_{famID}.tsv"
    params:
        info = "".join([f"--info-field {x} " for x in info_fields]),
        csq = "".join([f"--csq-column {x} " for x in csq_columns])
    container:
        "docker://brentp/slivar:v0.2.7"
    resources:
        nodes = 1
    shell:
        """
        slivar tsv \
            -s slivar_comphet \
               {params.info} \
               -i genic \
               -c CSQ \
               {params.csq} \
               -g slivar/lof_lookup.txt \
               -g slivar/clinvar_gene_desc.txt \
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
    resources:
        nodes = 1
    shell:
        """
        # get header from first file and drop it from other files
        # and make sure slivar_comphet id is unique
        #awk 'NR == FNR || FNR > 1 {{ sub(/^slivar_comphet/, "slivar_comphet_"NR, $0); print; }}' {input.filtered} {input.comphet}
        cat {input.filtered} {input.comphet} | sed '1 s/gene_description_1/lof/;s/gene_description_2/clinvar_trait/;' > {output}
        """
rule report:
    input: 
        rules.merge_tsv.output
    output:
        genic_tier1="reports/{famID}_genic_tier1.tsv",
        genic_tier2="reports/{famID}_genic_tier2.tsv",
        nongenic= "reports/{famID}_nongenic.tsv"
    script:
        "report_format.py"
