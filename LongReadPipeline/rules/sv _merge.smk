#!/usr/bin/env python
import glob
import os

configfile: "config.yaml"

# INPUT FASTQ folder
FASTQDIR = config["resources"]["fastqdir"]

max_distance = config["survivor_parameters"]["max_distance"]
min_call = config["survivor_parameters"]["min_call"]
svtype = config["survivor_parameters"]["svtype"]
strand = config["survivor_parameters"]["strand"]
min_size = config["survivor_parameters"]["min_size"]


######## RULES ##########  

rule samples:
    input:
        FASTQDIR
    output:
        sample_list = "{FASTQDIR}/sample_list.txt"
    shell:
         "printf '{FASTQDIR}/{sample}_lra/sv_calls/{sample}_lra_cutesv_filtered.vcf.gz\n{FASTQDIR}/{sample}_win/sv_calls/{sample}_win_cutesv_filtered.vcf.gz\n' > {output}"         


rule survivor:
    input:
        max_distance = max_distance
        min_call = min_call
        svtype = svtype
        strand = strand
        min_size = min_size
        sample_list = rules.samples.output.sample_list
    output:
        VCF = "{FASTQDIR}/{sample}_survivor/{sample}_cutesv_filtered.vcf"
    conda: "../envs/SV_env.yml"
    threads: thread_n
    shell:
        "SURVIVOR merge {input.sample_list} {input.max_distance} {input.min_call} {input.svtype} {input.strand} 0 {input.min_size} {output.VCF}"



rule sort_vcf:
    input:
        VCF = rules.survivor.output.VCF
    output:
        VCF = temp("{FASTQDIR}/{sample}/sv_calls/{sample}_cutesv_filtered.vcf")
    conda: "../envs/SV_env.yml"
    shell:
         "vcfsort {input.VCF} > {output.VCF}"



rule index_vcf:
    input:
         VCF = rules.sort_vcf.output.VCF
    output:
         VCF = "{FASTQDIR}/{sample}/sv_calls/{sample}_cutesv_filtered.vcf.gz"
    conda: "../envs/SV_env.yml"
    shell:
         "cat {input.VCF} | bgziptabix {output.VCF}"






