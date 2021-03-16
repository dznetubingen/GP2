#!/usr/bin/env python
import glob
import os

configfile: "config.yaml"

FASTQDIR = config["resources"]["fastqdir"]
VCF = config["resources"]["vcf"]
REF = config["resources"]["reference"]
BAM = "../pipeline-structural-variation/alignment/{sample}_lra.bam"
sample = config["resources"]["sample"]
rule phasing:
    input:
        bam = BAM,
        vcf = FASTQSHORT,
        ref = REF
    output:
        directory("{FASTQDIR}/whatshap_output")
    singularity:
        "../singularityIMG/long_read_tools.simg"
    threads: 15
    shell:
        """
	whatshap phase -o {output}/phased.vcf --reference={input.ref} {input.vcf} {input.bam}
	"""
