#!/bin/env python
import glob

configfile: "config.yaml"

FASTQDIR = config["resources"]["fastqdir"]

rule qc:
    input:
        FASTQDIR
    output:
        directory("{FASTQDIR}/Nanoplot")
    conda:
        "../envs/nanoplot.yml"
    threads: 30
    shell:
        """
        gunzip {input}/all_pass_reads.fastq.gz
        Nanoplot --fastq {input}/all_pass_reads.fastq -o {input}/Nanoplot
        """ 
