#!/bin/env python
import glob

configfile: "config.yaml"

FASTQDIR = config["resources"]["fastqdir"]


rule pychopper:
    input:
        FASTQDIR
    output:
        directory("{FASTQDIR}/pychopper_output")
    conda:
        "../envs/pychopper.yml"
    threads: 15
    shell:
        """
		cdna_classifier.py -t {threads} -x PCS109 -r {input}/pychopper_output/report.pdf -S {input}/pychopper_output/statistics.txt {input}/all_pass_reads.fastq {input}/pychopper_output/pychopper_classified_reoriented.fastq
		"""
