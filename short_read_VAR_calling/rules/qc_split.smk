#!/usr/bin/env python3
import glob
import os
import sys

configfile: "config.yaml"
FASTQDIR = config["resources"]["fastqdir"]
GENOMEDIR = config["resources"]["genomedir"]

## Parsing output directory
#OUTDIR = config["resources"]["outdir"]

SAMPLES, = glob_wildcards(FASTQDIR + "/{sample}_R1.fastq.gz")


singularity: "docker://continuumio/miniconda3"

localrules: multiqc

rule check_sample:
    run:
        print('These are sample names:')
        for sample in glob.glob(FASTQDIR + "/{sample}_R1.fastq.gz"):
            print(sample)

rule fastp_qc:
    input:
        R1 = FASTQDIR + "/{sample}_R1.fastq.gz",
        R2 = FASTQDIR + "/{sample}_R2.fastq.gz"
    output:
        R1_qc = "qc/{sample}_R1.fastq.gz",
        R2_qc = "qc/{sample}_R2.fastq.gz",
        JSON = "qc/{sample}_fastp.json"
    conda:
        "../envs/short_read_mapping.yml"
    params:
        qualified_quality_phred = 15,
        unqualified_percent_limit = 40,
        nova_seq = "-g"
    threads: 5
    log:
        "logs/fastp_{sample}.log"
    #resources:
    #     mem= "1gb",
    #     walltime= "04:00:00",
    #     node=1
    shell:
        """
        fastp --thread {threads} -i {input.R1} -o {output.R1_qc} -I {input.R2} -O {output.R2_qc} \
                -j {output.JSON} -q {params.qualified_quality_phred} \
                -u {params.unqualified_percent_limit} \
                {params.nova_seq} >/dev/null
        """

rule multiqc:
    input:
        expand("qc/{sample}_fastp.json", sample=SAMPLES)
    output:
        "qc/multiqc.html"
    conda:
        "../envs/short_read_mapping.yml"
    log: 
        "logs/multiqc.log"  
    shell:
        "multiqc -f qc -o qc 2> {log}" 

