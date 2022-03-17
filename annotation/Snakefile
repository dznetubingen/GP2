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
include: "rules/gene_panel.smk"


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
        "concordance/gp2_wgs_reduced_booster.pgen",
        expand("vcfs_per_fam/{famID}.tfam", famID=ped_dict.keys()),
        expand("reports/{famID}.xlsx", famID=ped_dict.keys()),
        expand("reports/{famID}_pubmed_table.xlsx", famID=ped_dict.keys()),
        expand("{panel}/{famID}.xlsx", panel= ['panel1','panel2'], famID=FAMIDS)

rule bcftools_norm:
    input:
        "inputs/gp2.vcf.gz"
    output:
        temp("vcfs/chrs/norm/{chromosome}.vcf.gz")
    params:
        chrom = "{chromosome}"
    threads:
        2
    conda:
        "envs/tools.yml"
    shell:
       """
       bcftools view -f PASS --regions {params.chrom} -Ou -s ^MH00724-TRO {input} \
       | bcftools norm \
        -f {config[ref_genome]} \
        --threads {threads} \
        -m - \
        -O u \
       | bcftools +fill-tags -O u -- -t AN,AC,FORMAT/VAF \
       | bcftools view --trim-alt-alleles --min-ac 1 -O v \
       | vcftools --vcf - \
                  --minGQ 20 \
                  --minDP 5 \
                  --recode \
                  --recode-INFO-all \
                  --stdout \
       | bcftools filter -S . -i 'VAF>0.2 && VAF < 0.75 | VAF < 0.02 | VAF > 0.98' -O z -o {output} \
        && tabix -p vcf --force {output}
       """


rule slivar_clean:
    input:
        vcf = rules.bcftools_norm.output
    output:
        "vcfs/chrs/filtered/{chromosome}.vcf.gz"
    container:
        "docker://zihhuafang/slivar_modified:0.2.7"
    threads:
        1
    resources:
        nodes = 1
    shell:
        """
        slivar expr --js /mnt/slivar/slivar-functions.js \
                    --pass-only \
                    --info 'variant.ALT[0] != "*" && variant.call_rate > 0.95' \
                    --vcf {input.vcf}  \
        | bgzip -c > {output} && tabix -p vcf --force {output}
        """

rule vep:
    input:
        vcf = rules.slivar_clean.output
    output:
        vcf = temp("vcfs/chrs/vep/{chromosome}.vcf.gz")
    container:
        "docker://zihhuafang/ensembl_vep_loftee:v105"
    threads: 10
    resources:
        nodes = 1
    shell:
        """
        vep \
        -i {input.vcf} \
        --ASSEMBLY GRCh38 \
        --buffer_size 50000 \
        --fasta {config[docker_ref_genome]} \
        --everything \
        --var_synonyms \
        --check_existing \
        --nearest symbol \
        --terms SO \
        --mirna \
        --check_existing \
        --clin_sig_allele 1 \
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
        --plugin TSSDistance \
        --plugin NMD \
        --plugin Downstream \
        --plugin DisGeNET,file=/opt/vep/.vep/vep_annotation/all_variant_disease_pmid_associations_final.tsv.gz,disease=1 \
        --plugin dbNSFP,/opt/vep/.vep/vep_annotation/dbNSFP4.2a_grch38.gz,Ensembl_transcriptid,MetaRNN_score,LRT_score,GERP++_RS,FATHMM_score,fathmm-MKL_coding_score,DANN_score,REVEL_score,PrimateAI_score,Geuvadis_eQTL_target_gene,gnomAD_genomes_controls_and_biobanks_nhomalt \
        --plugin LoF,loftee_path:/opt/vep/.vep/Plugins/loftee_hg38,\
human_ancestor_fa:/opt/vep/.vep/Plugins/loftee_hg38/human_ancestor.fa.gz,\
conservation_file:/opt/vep/.vep/Plugins/loftee_hg38/loftee.sql,\
gerp_bigwig:/opt/vep/.vep/Plugins/loftee_hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
        --custom {config[clinvar]},ClinVar,vcf,exact,0,CLNSIG,CLNDN \
        -o {output.vcf} \
        --compress_output bgzip \
        --vcf \
        --stats_text
        """

#rule split_vep:
#    input:
#       vcf = rules.vep.output
#    output:
#        temp("vcfs/chrs/split_vep/{chromosome}.vcf.gz")
#    conda:
#        "envs/tools.yml"
#    threads: 1
#    resources:
#        nodes = 1
#    shell:
#        """
#        bcftools +split-vep \
#          -p vep \
#          -a CSQ \
#          -c MaxEntScan_alt:Float,MaxEntScan_diff:Float,MaxEntScan_ref:Float,SpliceAI_pred_DP_AG:Float,SpliceAI_pred_DP_AL:Float,SpliceAI_pred_DP_DG:Float,SpliceAI_pred_DP_DL:Float,SpliceAI_pred_DS_AG:Float,SpliceAI_pred_DS_AL:Float,SpliceAI_pred_DS_DG:Float,SpliceAI_pred_DS_DL:Float,SpliceAI_pred_SYMBOL:String,CLIN_SIG:String \
#          {input.vcf} | \
#        sed -e 's/vepMaxEntScan_alt,Number=./vepMaxEntScan_alt,Number=1/g' \
#        -e 's/vepMaxEntScan_diff,Number=./vepMaxEntScan_diff,Number=1/g'  \
#        -e 's/vepMaxEntScan_ref,Number=./vepMaxEntScan_ref,Number=1/g'   \
#        -e 's/vepSpliceAI_pred_DP_AG,Number=./vepSpliceAI_pred_DP_AG,Number=1/g'  \
#        -e 's/vepSpliceAI_pred_DP_AL,Number=./vepSpliceAI_pred_DP_AL,Number=1/g'  \
#        -e 's/vepSpliceAI_pred_DP_DG,Number=./vepSpliceAI_pred_DP_DG,Number=1/g'  \
#        -e 's/vepSpliceAI_pred_DP_DL,Number=./vepSpliceAI_pred_DP_DL,Number=1/g'  \
#        -e 's/vepSpliceAI_pred_DS_AG,Number=./vepSpliceAI_pred_DS_AG,Number=1/g'  \
#        -e 's/vepSpliceAI_pred_DS_AL,Number=./vepSpliceAI_pred_DS_AL,Number=1/g'  \
#        -e 's/vepSpliceAI_pred_DS_DG,Number=./vepSpliceAI_pred_DS_DG,Number=1/g'  \
#        -e 's/vepSpliceAI_pred_DS_DL,Number=./vepSpliceAI_pred_DS_DL,Number=1/g' \
#         | bgzip -c > {output} && tabix -p vcf --force {output}
#        """

rule prep_vcf_CADD:
    input:
        rules.slivar_clean.output
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
        vcf = rules.slivar_clean.output,
        tfam = "inputs/gp2_all.tfam"
    output:
        temp("vcfs/gcount/{chromosome}.vcf")
    conda:
        "envs/tools.yml"
    threads: 1
    resources:
        nodes = 1
    shell:
        """
        #get n of homalt, n of het and allele count in cases and controls
        SnpSift caseControl -tfam {input.tfam} {input.vcf} > {output}

        """


rule add_anno_vcf:
    input:
        vcf = rules.vep.output,
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
        slivar expr --js /mnt/slivar/slivar-functions.js \
                    --pass-only \
                    -g /mnt/slivar/TOPMed_freeze8_PASS.zip \
                    -g /mnt/slivar/gnomad_v3.1.2.zip \
                    --info 'variant.ALT[0] != "*" && variant.call_rate > 0.95' \
                    --vcf {input.vcf}  \
        | bgzip -c > {output} && tabix -p vcf --force {output}
        """

rule plink_bfiles:
    input:
        vcf= rules.slivar_clean.output
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
             --keep-allele-order \
             --set-all-var-ids @:#:\$r:\$a \
             --double-id \
             --new-id-max-allele-len 300 missing \
             --make-bed \
             --out bfiles/{params.chr}
        """

rule Gathervcf:
    input:
        lambda wildcards: expand("final_vcfs/{chromosome}.vcf.gz", chromosome=chr_id)
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
        vcf = "final_vcfs/GP2_annotate.vcf.gz"
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
             --set-all-var-ids @:#:\$r:\$a \
             --new-id-max-allele-len 300 missing \
             --keep-allele-order \
             --double-id \
             --make-bed \
             --out bfiles/all_chrs_merged
        """


rule plink_pfiles:
    input:
        vcf= "final_vcfs/GP2_annotate.vcf.gz"
    output:
        multiext("concordance/all_chrs_merged", ".pgen", ".psam", ".pvar")
    conda:
        "envs/tools.yml"
    threads: 5
    resources:
        nodes = 1
    shell:
        """
        plink2 --vcf {input.vcf} \
             --threads {threads} \
             --set-all-var-ids @:#:\$r:\$a \
             --new-id-max-allele-len 300 missing \
             --double-id \
             --make-pgen \
             --out concordance/all_chrs_merged
        """

rule plink2_wgs_reduced:
    input:
        pgen= "concordance/all_chrs_merged.pgen",
        booster_snp= "concordance/extract_booster_snp.txt"
    output:
        multiext("concordance/gp2_wgs_reduced_booster", ".pgen", ".psam", ".pvar")
    conda:
        "envs/tools.yml"
    threads: 5
    params: "concordance/all_chrs_merged"
    resources:
        nodes = 1
    shell:
        """
        plink2 --pfile {params} \
             --threads {threads} \
             --extract {input.booster_snp} \
             --double-id \
             --make-pgen \
             --out concordance/gp2_wgs_reduced_booster
        """
