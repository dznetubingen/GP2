#!/usr/bin/env python
import glob
import os

configfile: "config.yaml"
FAST5DIR=config["resources"]["fast5dir"]
FASTQDIR = config["resources"]["fastqdir"]


localrules: zcat

rule guppy:
    input:
        FAST5DIR
    output:
        FASTQDIR
    params:
        config_file= config["resources"]["guppy_config"]
    threads:
        15
    resources:
        nvidia_gpu=1
    singularity:
        "../singularityIMG/guppy-gpu-4.4.1.simg"
    shell:
        """
        guppy_basecaller --compress_fastq -i {input} --recursive \
                --save_path {output} \
                --num_callers {threads} \
                -c {params.config_file} \
                --trim_strategy none \
                --device auto \
                --qscore_filtering --min_qscore 7
        """
        
            
rule zcat:
    input:
        FASTQDIR
    output:
        "{FASTQDIR}/Q7_pass_reads.fastq.gz"
    threads:
        5
    shell:
        """
        zcat {input}/pass/*.fastq.gz > {output}
        """
