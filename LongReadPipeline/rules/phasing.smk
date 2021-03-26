#!/usr/bin/env python
import glob
import os

configfile: "config.yaml"

FASTQDIR = config["resources"]["fastqdir"]
VCF = config["resources"]["vcf"]
REF = config["resources"]["reference"]
BAM = "{FASTQDIR}/{sample}_win/alignment/{sample}_win.bam"
sample = config["resources"]["sample"]
rule phasing:
    input:
        bam = BAM,
        vcf = FASTQSHORT,
        ref = REF
    output:
        directory("{FASTQDIR}/whatshap_output")
    singularity:
        "docker://avalillarionova/lrp_tools:latest"
    threads: 15
    shell:
        """
	whatshap phase -o {output}/phased_winnowmap.vcf --reference={input.ref} {input.vcf} {input.bam}
	"""
