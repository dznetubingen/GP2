#!/usr/bin/env python
import glob
import os

configfile: "config.yaml"

FAST5DIR=config["resources"]["fast5dir"]
FASTQDIR = config["resources"]["fastqdir"]
BAM = "{FASTQDIR}/{sample}_win/alignment/{sample}_win.bam"
REF = config["resources"]["reference"]

rule phasing:
    input:
        fastq = FASTQDIR,
        fast5 = FAST5DIR,
        bam = BAM,
        ref = REF
    output:
        directory("{FASTQDIR}/nanopolish_output")
    singularity:
        "docker://avalillarionova/lrp_tools:latest"
    threads: 15
    shell:
        """
	nanopolish index -d {input.fast5}/ {input.fastq}/Q7_pass_reads.fastq.gz
	nanopolish call-methylation -t {threads} -r {input.fastq}/Q7_pass_reads.fastq.gz -b {input.bam} -g {input.ref} > {output}/methylation_calls.tsv
	"""
