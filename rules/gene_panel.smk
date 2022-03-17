#!/usr/bin/env python
import glob
import pandas as pd
from pathlib import Path
from snakemake.utils import min_version

shell.prefix("set -o pipefail; umask 002; ")  # set g+w

if not os.path.exists("panel"):
    os.makedirs("panel")

#Snakemake config
min_version("5.5")
configfile: "configs/config.yaml"

#parse wildcards
df= pd.read_csv(config["gp2_ped"], sep='\t', header=None, names=['family_id','id','parental_id','maternal_id','sex','phenotype'])
FAMIDS = list(df.loc[df.phenotype==2, 'family_id'].unique())

gene_descs= [
    '/mnt/slivar/pLI_lookup.txt',
    '/mnt/slivar/oe_lof_upper_lookup.txt',
    '/mnt/slivar/oe_mis_upper_lookup.txt',
    '/mnt/slivar/oe_syn_upper_lookup.txt',
    '/mnt/slivar/clinvar_gene_desc.txt',
    '/mnt/slivar/hgncSymbol.inheritance.tsv',
    '/mnt/slivar/lookup_description_symbol.txt'
    ]


rule panel1:
    input:
        vcf = "final_vcfs/GP2_annotate.vcf.gz",
        ped = "inputs/gp2_all.ped", 
        bed = "inputs/PD_panel_all_hg38.bed"   #include both core and extensive panels
    output:
        "panel1/GP2.vcf"
    container:
        "docker://zihhuafang/slivar_modified:0.2.7"
    params:
        gene_desc = "".join([f"-g {x} " for x in gene_descs])
    threads: 1
    resources:
        nodes = 1
    shell:
        """
        slivar expr --vcf {input.vcf} \
                    -o {output} \
                    --ped {input.ped} \
                    --js slivar/slivar-functions.js \
                    --pass-only \
                    --info 'INFO.gnomad_AF < 0.05 && variant.ALT[0] != "*" && variant.call_rate > 0.95' \
                    --region {input.bed} \
                    --family-expr "HOM:fam.every(function(s) {{return s.hom_alt == s.affected && hq1(s)}})" \
                    --family-expr "HET:fam.every(function(s) {{return s.het == s.affected && hq1(s) }})"
        """

rule panel2:
    input:
        vcf = "final_vcfs/GP2_annotate.vcf.gz",
        ped = "inputs/gp2_all.ped",
        bed = "inputs/parkinson_panel_disease_gp2.bed"   #panel proposed by Ignacio
    output:
        "panel2/GP2.vcf"
    container:
        "docker://zihhuafang/slivar_modified:0.2.7"
    params:
        gene_desc = "".join([f"-g {x} " for x in gene_descs])
    threads: 1
    resources:
        nodes = 1
    shell:
        """
        slivar expr --vcf {input.vcf} \
                    -o {output} \
                    --ped {input.ped} \
                    --js slivar/slivar-functions.js \
                    --pass-only \
                    --info 'INFO.gnomad_AF < 0.05 && variant.ALT[0] != "*" && variant.call_rate > 0.95' \
                    --region {input.bed} \
                    --family-expr "HOM:fam.every(function(s) {{return s.hom_alt == s.affected && hq1(s)}})" \
                    --family-expr "HET:fam.every(function(s) {{return s.het == s.affected && hq1(s) }})"
        """

info_fields = [
    'gnomad_AF',
    'gnomad_popmax_af',
    'gnomad_nhomalt',
    'gnomad_AC',
    'TOPMed8_AF',
    'CADD_PHRED',
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
    'miRNA',
    'NMD',
    'DisGeNET'
    ]

rule panel_report:
    input:
       vcf= "{panel}/GP2.vcf",
       ped = "inputs/gp2_all.ped"
    output:
        "{panel}/all_samples.tsv"
    container:
        "docker://zihhuafang/slivar_modified:0.2.7"
    params:
        info = "".join([f"--info-field {x} " for x in info_fields]),
        csq = "".join([f"--csq-column {x} " for x in csq_columns]),
        gene_desc = "".join([f"-g {x} " for x in gene_descs])
    threads: 1
    resources:
        nodes = 1
    shell:
        """
        slivar tsv \
               {params.info} \
               -s HOM \
               -s HET \
               -c CSQ \
               {params.csq} \
               {params.gene_desc} \
               -p {input.ped} \
               {input.vcf} > {output} 
        """

rule fix_tsv:
    input:
        "{panel}/all_samples.tsv"
    output:
        "{panel}/fixed_all_samples.tsv"
    threads: 1
    resources:
        nodes = 1
    shell:
        """
        cat {input} \
        | sed '1 s/gene_description_1/gnomAD_pLI/;s/gene_description_2/gnomAD_oe_lof_CI90/;s/gene_description_3/gnomAD_oe_mis_CI90/;s/gene_description_4/gnomAD_oe_syn_CI90/;s/gene_description_5/clinvar_gene_description/;s/gene_description_6/MOI/;s/gene_description_7/gene_fullname/;' > {output}
        """

rule format_report:
    input:
        tsv = "panel1/fixed_all_samples.tsv",
        core_panel = "inputs/PD_core_panel.list",
        extend_panel = "inputs/PD_extended_panel.list"
    output:
        expand("panel1/{famID}.xlsx", famID=FAMIDS)
    threads: 1
    script:
        ### modify the report script
        "{workflow.basedir}/scripts/report_panel.py"


rule format_report2:
    input:
        tsv = "panel2/fixed_all_samples.tsv"
    output:
        expand("panel2/{famID}.xlsx", famID=FAMIDS)
    threads: 1
    script:
        ### modify the report script
        "{workflow.basedir}/scripts/report_panel2.py"
