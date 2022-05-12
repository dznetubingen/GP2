#!/usr/bin/env python
import glob
import os
from snakemake.utils import min_version


#log output dir
if not os.path.exists("slurm"):
    os.makedirs("slurm")

shell.prefix("set -o pipefail; umask 002; ")  # set g+w

#Snakemake config
min_version("5.5")
configfile: "configs/config.yaml"

rule plink_bfile:
    input:
        vcf = "final_vcfs/GP2_annotate.vcf.gz"
    output:
        "bfiles/all_chrs_merged.bed"
    conda:
        "../envs/tools.yml"
    threads: 5
    resources:
        nodes = 1
    shell:
        """
        plink2 --vcf {input.vcf} \
             --threads {threads} \
             --set-missing-var-ids @:#:\$r:\$a \
             --new-id-max-allele-len 500 missing \
             --keep-allele-order \
             --double-id \
             --make-bed \
             --out bfiles/all_chrs_merged
        """

rule plink_pfile:
    input:
        vcf= "final_vcfs/GP2_annotate.vcf.gz"
    output:
        multiext("concordance/all_chrs_merged", ".pgen", ".psam", ".pvar")
    conda:
        "../envs/tools.yml"
    threads: 5
    resources:
        nodes = 1
    shell:
        """
        plink2 --vcf {input.vcf} \
             --threads {threads} \
             --double-id \
             --make-pgen \
             --out concordance/all_chrs_merged
        """

# need to remake "concordance/extract_booster_snp.txt" with hg38 vcf
rule plink2_wgs_reduced:
    input:
        pgen= "concordance/all_chrs_merged.pgen",
        booster_snp= "concordance/extract_booster_snp.txt"
    output:
        multiext("concordance/gp2_wgs_reduced_booster", ".pgen", ".psam", ".pvar")
    conda:
        "../envs/tools.yml"
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

rule merge_with_1kg:
    input:
        vcf= "final_vcfs/GP2_annotate.vcf.gz",
        ref_panel = "/mnt/vol009/fangz/ref_panel/gp2_panel_hg38/gp2_ref_panel_hg38.vcf.gz"
    output:
        "qc/gp2_merged_refpanel_to_qc.vcf.gz"
    conda:
        "../envs/tools.yml"
    threads: 5
    resources:
        nodes = 1
    shell:
        """
        bcftools merge -m none -O v {input.ref_panel} {input.vcf}\
        | vcftools --vcf - \
           --maf 0.05 \
           --hwe 0.0001 \
           --max-missing 0.95 \
           --recode \
           --stdout \
        | bgzip -c > {output} 

        tabix -p vcf -f {output}
        """

rule prep_bed_qc:
    input:
       rules.merge_with_1kg.output
    output:
        multiext("qc/gp2_merged_refpanel_to_qc", ".bed", ".bim", ".fam")
    conda:
        "../envs/tools.yml"
    threads: 5
    params: "qc/gp2_merged_refpanel_to_qc"
    resources:
        nodes = 1
    shell:
        """
        plink2 --vcf {input} \
             --threads {threads} \
             --double-id \
             --make-bed \
             --out {params}
        """
