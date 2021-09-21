#!/usr/bin/env python
import glob
import os
from snakemake.utils import min_version


#log output dir
if not os.path.exists("slurm"):
    os.makedirs("slurm")

shell.prefix("set -o pipefail; umask 002; ")  # set g+w

##include rules
include: "rules/CADD.smk"
include: "rules/slivar_filter.smk"

#Snakemake config
min_version("5.5")
configfile: "configs/config.yaml"

chr_list = list(range(1,23))+["X", "Y"]
chr_id= ["chr" + str(i) for i in chr_list]

rule all:
    input:
        expand("bfiles/{chromosome}.bed", chromosome = chr_id),
        expand("CADD/{chromosome}.tsv.gz", chromosome = chr_id),
        expand("final_vcfs/{chromosome}.vcf.gz", chromosome = chr_id),
        "final_vcfs/GP2_annotate.vcf.gz",
        "bfiles/all_chrs_merged.bed",
        "short_list/all_chrs.tsv"

rule bcftools_norm:
    input:
        "gp2.vcf.gz"
    output:
        temp("gp2_norm.vcf.gz")
    threads:
        20
    conda:
        "envs/tools.yml"
    shell:
       """
       bcftools norm \
        -f {config[ref_genome]} \
        --threads {threads} \
        -m - \
        -O z \
        --output {output} \
        {input} && tabix -p vcf --force {output} 
       """

rule split_by_chr:
    input:
        vcf= rules.bcftools_norm.output
    output:
        temp("vcfs_chr/{chromosome}.vcf.gz")
    params:
        chr = "{chromosome}"
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
    threads: 1
    resources:
        nodes = 1
    shell:
        """
        gatk SelectVariants \
          -R {config[ref_genome]} \
          -V {input.vcf} \
          -L {params.chr} \
          --exclude-filtered true \
          -O {output}
        """

rule vep:
    input:
        vcf = "vcfs_chr/{chromosome}.vcf.gz"
    output:
        vcf = temp("vcf_vep/{chromosome}.vcf.gz")
    container:
        "docker://ensemblorg/ensembl-vep:release_104.3"
    threads: 10
    resources:
        nodes = 1
    shell:
        """
        vep \
        -i {input.vcf} \
        --ASSEMBLY GRCh38 \
        --buffer_size 50000 \
        --polyphen b \
        --fasta {config[ref_genome]} \
        --sift b \
        --biotype \
        --hgvs \
        --protein \
        --domains \
        --pubmed \
        --ccds \
        --check_existing \
        --nearest symbol \
        --gene_phenotype \
        --canonical \
        --regulatory \
        --terms SO \
        --appris \
        --check_existing \
        --clin_sig_allele 1 \
        --var_synonyms \
        --variant_class \
        --af --max_af --af_1kg --af_esp --af_gnomad \
        --force_overwrite \
        --cache \
        --pick \
        --per_gene \
        --offline \
        --dir_cache /opt/vep/.vep  \
        --dir_plugins /opt/vep/.vep/Plugins/ \
        --plugin UTRannotator,./vep_annotation/uORF_5UTR_GRCh38_PUBLIC.txt \
        --plugin SpliceAI,snv=./vep_annotation/spliceai_scores.raw.snv.hg38.vcf.gz,indel=./vep_annotation/spliceai_scores.raw.indel.hg38.vcf.gz \
        --plugin SpliceRegion,Extended \
        --plugin MaxEntScan,/opt/vep/.vep/Plugins/maxEntScan \
        --plugin LoFtool,./vep_annotation/LoFtool_scores.txt \
        --plugin TSSDistance \
        --plugin DisGeNET,file=./vep_annotation/all_variant_disease_pmid_associations_final.tsv.gz \
        --plugin dbNSFP,./vep_annotation/dbNSFP4.2a_grch38.gz,MetaRNN_score,LRT_score,GERP++_RS,FATHMM_score,fathmm-MKL_coding_score,DANN_score,REVEL_score,PrimateAI_score,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_MedGen_id,clinvar_OMIM_id,Geuvadis_eQTL_target_gene,gnomAD_genomes_controls_and_biobanks_nhomalt \
        -o {output.vcf} \
        --compress_output bgzip \
        --fork {threads} \
        --vcf \
        --stats_text
        """

rule split_vep:
    input:
       vcf = rules.vep.output
    output:
        temp("split_vep/{chromosome}.vcf.gz")
    conda:
        "envs/tools.yml"
    threads: 5
    resources:
        nodes = 1
    shell:
        """
        bcftools +split-vep \
          -p vep \
          -a CSQ \
          -c MaxEntScan_alt:Float,MaxEntScan_diff:Float,MaxEntScan_ref:Float,SpliceAI_pred_DP_AG:Float,SpliceAI_pred_DP_AL:Float,SpliceAI_pred_DP_DG:Float,SpliceAI_pred_DP_DL:Float,SpliceAI_pred_DS_AG:Float,SpliceAI_pred_DS_AL:Float,SpliceAI_pred_DS_DG:Float,SpliceAI_pred_DS_DL:Float,SpliceAI_pred_SYMBOL:String,CLIN_SIG:String \
          {input.vcf} | \
        sed -e 's/vepMaxEntScan_alt,Number=./vepMaxEntScan_alt,Number=1/g' \
        -e 's/vepMaxEntScan_diff,Number=./vepMaxEntScan_diff,Number=1/g'  \
        -e 's/vepMaxEntScan_ref,Number=./vepMaxEntScan_ref,Number=1/g'   \
        -e 's/vepSpliceAI_pred_DP_AG,Number=./vepSpliceAI_pred_DP_AG,Number=1/g'  \
        -e 's/vepSpliceAI_pred_DP_AL,Number=./vepSpliceAI_pred_DP_AL,Number=1/g'  \
        -e 's/vepSpliceAI_pred_DP_DG,Number=./vepSpliceAI_pred_DP_DG,Number=1/g'  \
        -e 's/vepSpliceAI_pred_DP_DL,Number=./vepSpliceAI_pred_DP_DL,Number=1/g'  \
        -e 's/vepSpliceAI_pred_DS_AG,Number=./vepSpliceAI_pred_DS_AG,Number=1/g'  \
        -e 's/vepSpliceAI_pred_DS_AL,Number=./vepSpliceAI_pred_DS_AL,Number=1/g'  \
        -e 's/vepSpliceAI_pred_DS_DG,Number=./vepSpliceAI_pred_DS_DG,Number=1/g'  \
        -e 's/vepSpliceAI_pred_DS_DL,Number=./vepSpliceAI_pred_DS_DL,Number=1/g' \
         | bgzip -c > {output} && tabix -p vcf --force {output}
        """

rule prep_vcf_CADD:
    input:
        vcf = "vcfs_chr/{chromosome}.vcf.gz",
    output:
        temp("CADD_input/{chromosome}.vcf")
    threads: 1
    resources:
        nodes = 1
    shell:
        """
        zcat {input} | awk '{{gsub(/^chr/,""); print}}' > {output}
        """

rule geno_count:
    input:
        vcf= "vcfs_chr/{chromosome}.vcf.gz"
    output:
        "gp2_gcount/{chromosome}.gcount.gz"
    params:
        chr = "{chromosome}"
    conda:
        "envs/tools.yml"
    threads: 5
    resources:
        nodes = 1
    shell:
        """
        plink2 --vcf {input.vcf} \
             --threads {threads} \
             --geno-counts cols=chrom,pos,ref,alt,homref,refalt,altxy,missing \
             --out gp2_gcount/{params.chr} &&
         bgzip gp2_gcount/{params.chr}.gcount && tabix -s1 -b2 -e2 --force {output}     
        """


rule add_anno_vcf:
    input:
        vcf = "split_vep/{chromosome}.vcf.gz",
        vcfanno_conf = "vcfanno_conf/vcfanno_{chromosome}.conf",
        CADD_score= "CADD/{chromosome}.tsv.gz",
        geno_count="gp2_gcount/{chromosome}.gcount.gz"
    output:
        "vcfs_vcfanno/{chromosome}.vcf.gz"
    conda:
        "envs/tools.yml"
    threads: 5
    resources:
        nodes = 1
    shell:
        """
        vcfanno -p {threads} {input.vcfanno_conf} {input.vcf} \
        | bgzip -c > {output} && tabix -p vcf --force {output}
        """

rule slivar_gnotate:
    input:
        vcf = "vcfs_vcfanno/{chromosome}.vcf.gz"
    output:
        temp("final_vcfs/{chromosome}.vcf")
    container:
        "docker://brentp/slivar:v0.2.5"
    threads: 1
    resources:
        nodes = 1
    shell:
        """ 
        slivar expr --js slivar/slivar-functions.js \
                    --pass-only \
                    -g slivar/gnomad.hg38.genomes.v3.fix.zip \
                    -g slivar/topmed.hg38.dbsnp.151.zip \
                    --info 'variant.ALT[0] != "*"' \
                    --vcf {input.vcf}  \
                    -o  {output}
        """

rule tabix: 
    input: 
        "final_vcfs/{chromosome}.vcf"
    output:
        "final_vcfs/{chromosome}.vcf.gz"
    conda:
        "envs/tools.yml"
    threads: 1
    resources:
        nodes = 1
    shell:
        """
        bgzip {input} && tabix -p vcf {output}
        """

rule plink_bfiles:
    input:
        vcf= "vcfs_chr/{chromosome}.vcf.gz"
    output:
        "bfiles/{chromosome}.bed"
    params:
        chr = "{chromosome}"
    conda:
        "envs/tools.yml"
    threads: 5
    resources:
        nodes = 1
    shell:
        """
        plink2 --vcf {input.vcf} \
             --threads {threads} \
             --make-bed \
             --out bfiles/{params.chr}
        """

rule Gathervcf:
    input:
        expand("final_vcfs/{chromosome}.vcf.gz", chromosome=chr_id)
    output:
        "final_vcfs/GP2_annotate.vcf.gz"
    conda:
        "envs/tools.yml"
    threads: 10
    resources:
        nodes = 1
    shell:
        """
        bcftools concat {input} \
                 --threads {threads} \
                 -O z \
                 -o {output} && tabix -p vcf {output}
         
        """

rule plink_gather:
    input:
        vcf= "final_vcfs/GP2_annotate.vcf.gz"
    output:
        "bfiles/all_chrs_merged.bed"
    conda:
        "envs/tools.yml"
    threads: 5
    resources:
        nodes = 1
    shell:
        """
        plink2 --vcf {input.vcf} \
             --threads {threads} \
             --make-bed \
             --out bfiles/all_chrs_merged
        """
