#!/usr/bin/env python
import glob
import os

configfile: "config.yaml"

FASTQDIR = config["resources"]["fastqdir"]
#RM = config["resources"]["reads_manifest"]
REF = config["resources"]["reference"]
ANNOT = config["resources"]["annotation"]
rule pychopper:
    input:
        fastqdir = FASTQDIR,
        readm = RM,
        ref = REF,
        annot = ANNOT
    output:
        directory("{FASTQDIR}/flair_output")
    conda:
        "../envs/flair.yml"
    threads: 15
    shell:
        """
	python flair/flair.py 123 -r {input.fastqdir}/pychopper_output/pychopper_output.fq -g {input.ref} -f {input.annot} --temp_dir {input.fastqdir}/temp_flair -o flair
	"""
