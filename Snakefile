#!/usr/bin/env python
import glob
import os

#include: "rules/qc_split.smk"
#include: "rules/align.smk"
include: "rules/haplotype_caller.smk"

#create tmpdir by yourself as $TMPDIR is shared with all the jobs on the drive, sometimes there's not enough space
if not os.path.exists("tmp_scratch"):
    os.makedirs("tmp_scratch")

#create log folder
if not os.path.exists("logs"):
    os.makedirs("logs")


configfile: "config.yaml"
FASTQDIR = config["resources"]["fastqdir"]
GENOMEDIR = config["resources"]["genomedir"]

SAMPLES, = glob_wildcards(FASTQDIR + "/{sample}_R1.fastq.gz")
chr_list = list(range(1,22))+["X", "Y", "M"]
chr_id= ["chr" + str(i) for i in chr_list]


singularity: "docker://continuumio/miniconda3"
localrules: multiqc, all

rule all:
    input:
#        "qc/multiqc.html",
#        expand("alignment/{sample}.bam", sample=SAMPLES),
#        expand("stats/{sample}_aln.stats", sample=SAMPLES),
#        expand("stats/{sample}.metrics.txt", sample=SAMPLES),
#        expand("bam_recal/{sample}_recal.table", sample=SAMPLES)
        expand("raw_vcf/{chromosome}.vcf.gz",chromosome=chr_id, sample=SAMPLES)
        
