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

rule plink_gather:
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
        "../envs/tools.yml"
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
