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
include: "rules/filter_per_fam.smk"

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
#        "concordance/gp2_wgs_reduced_booster.pgen"
        expand("reports/{famID}_{varlist}.tsv", famID=FAM, varlist=['genic_tier1','genic_tier2','nongenic'])

rule bcftools_norm:
    input:
        "inputs/gp2.vcf.gz"
    output:
        temp("vcfs/gp2_norm.vcf.gz")
    threads:
        20
    conda:
        "envs/tools.yml"
    shell:
       """
       bcftools view -f PASS -Ou -s ^MH00724-TRO {input} \
       | bcftools norm \
        -f {config[ref_genome]} \
        --threads {threads} \
        -m - \
        -O z \
        --output {output} \
        && tabix -p vcf --force {output} 
       """

rule split_by_chr:
    input:
        vcf= rules.bcftools_norm.output
    output:
        temp("vcfs/per_chr/{chromosome}.vcf.gz")
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
          --exclude-non-variants true \
          -O {output}
        """

rule vep:
    input:
        vcf = rules.split_by_chr.output
    output:
        vcf = temp("vcfs/per_chr/vep/{chromosome}.vcf.gz")
    container:
        "docker://ensemblorg/ensembl-vep:release_105.0"
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
        --plugin UTRannotator,/opt/vep/.vep/vep_annotation/uORF_5UTR_GRCh38_PUBLIC.txt \
        --plugin SpliceAI,snv=/opt/vep/.vep/vep_annotation/spliceai_scores.masked.indel.hg38.vcf.gz,indel=/opt/vep/.vep/vep_annotation/spliceai_scores.masked.indel.hg38.vcf.gz \
        --plugin SpliceRegion,Extended \
        --plugin MaxEntScan,/opt/vep/.vep/Plugins/maxEntScan \
        --plugin LoFtool,/opt/vep/.vep/vep_annotation/LoFtool_scores.txt \
        --plugin TSSDistance \
        --plugin NMD \
        --plugin DisGeNET,file=/opt/vep/.vep/vep_annotation/all_variant_disease_pmid_associations_final.tsv.gz \
        --plugin dbNSFP,/opt/vep/.vep/vep_annotation/dbNSFP4.2a_grch38.gz,Ensembl_transcriptid,MetaRNN_score,LRT_score,GERP++_RS,FATHMM_score,fathmm-MKL_coding_score,DANN_score,REVEL_score,PrimateAI_score,Geuvadis_eQTL_target_gene,gnomAD_genomes_controls_and_biobanks_nhomalt \
        --custom {config[clinvar]},ClinVar,vcf,exact,0,CLNSIG,CLNDN \
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
        temp("vcfs/per_chr/split_vep/{chromosome}.vcf.gz")
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
        vcf = rules.split_by_chr.output,
    output:
        temp("CADD_input/{chromosome}.vcf")
    threads: 1
    resources:
        nodes = 1
    shell:
        """
        zcat {input} | awk '{{gsub(/^chr/,""); print}}' > {output}
        """

rule gcount_case_control:
    input:
        vcf = rules.split_by_chr.output,
        tfam = "input_vcfs/gp2_all.tfam"
    output:
        temp("vcfs/gcount/{chromosome}.vcf")
    conda:
        "envs/tools.yml"
    resources:
        nodes = 1
    shell:
        """
        #get n of homalt, n of het and allele count in cases and controls
        SnpSift caseControl -tfam {input.tfam} {input.vcf} > {output}

        """


rule add_anno_vcf:
    input:
        vcf = rules.split_vep.output,
        vcfanno_conf = "vcfanno_conf/vcfanno_{chromosome}.conf",
        CADD_score= "CADD/{chromosome}.tsv.gz",
        gcount = rules.gcount_case_control.output
    output:
        count_vcf= temp("vcfs/gcount/{chromosome}.vcf.gz"),
        anno_vcf=temp("vcfs/vcfanno/{chromosome}.vcf.gz")
    conda:
        "envs/tools.yml"
    threads: 5
    resources:
        nodes = 1
    shell:
        """ 
        bgzip {input.gcount} && tabix -p vcf --force {output.count_vcf} &&
        vcfanno -p {threads} {input.vcfanno_conf} {input.vcf} \
        | bgzip -c > {output.anno_vcf} && tabix -p vcf --force {output.anno_vcf}
        """

rule slivar_gnotate:
    input:
        vcf = "vcfs/vcfanno/{chromosome}.vcf.gz"
    output:
        "final_vcfs/{chromosome}.vcf.gz"
    container:
        "docker://zihhuafang/slivar_modified:0.2.7"
    threads: 
        1 
    resources:
        nodes = 1
    shell:
        """ 
        slivar expr --js slivar/slivar-functions.js \
                    --pass-only \
                    -g slivar/gnomad.hg38.v3.custom.zip \
                    -g slivar/topmed.hg38.dbsnp.151.zip \
                    --info 'variant.ALT[0] != "*" && variant.call_rate > 0.95' \
                    --vcf {input.vcf}  \
        | bgzip -c > {output} && tabix -p vcf --force {output}
        """


rule plink_bfiles:
    input:
        vcf= rules.split_by_chr.output
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
             --set-all-var-ids @:#:\$r-\$a \
             --new-id-max-allele-len 50 missing \
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
                 -o {output} && tabix -p vcf --force {output}
         
        """

rule plink_gather:
    input:
        vcf= "final_vcfs/GP2_annotate.vcf.gz"
    output:
        "bfiles/all_chrs_merged.bed"
    conda:
        "envs/tools.yml"
    threads: 20
    resources:
        nodes = 1
    shell:
        """
        plink2 --vcf {input.vcf} \
             --threads {threads} \
             --set-all-var-ids @:#:\$r-\$a \
             --new-id-max-allele-len 50 missing \
             --make-bed \
             --out bfiles/all_chrs_merged
        """


rule plink_pfiles:
    input:
        vcf= "final_vcfs/GP2_annotate.vcf.gz"
    output:
        "concordance/all_chrs_merged.pgen"
    conda:
        "envs/tools.yml"
    threads: 20
    resources:
        nodes = 1
    shell:
        """
        plink2 --vcf {input.vcf} \
             --threads {threads} \
             --set-all-var-ids @:#-\$r-\$a \
             --new-id-max-allele-len 50 missing \
             --make-pgen \
             --out concordance/all_chrs_merged
        """

rule plink2_wgs_reduced:
    input:
        pgen= "concordance/all_chrs_merged.pgen",
        booster_snp= "concordance/extract_booster_snp.txt"
    output:
        "concordance/gp2_wgs_reduced_booster.pgen"
    conda:
        "envs/tools.yml"
    threads: 20
    params: "concordance/all_chrs_merged"
    resources:
        nodes = 1
    shell:
        """
        plink2 --pfile {params} \
             --threads {threads} \
             --extract {input.booster_snp} \
             --make-pgen \
             --out concordance/gp2_wgs_reduced2booster
        """
