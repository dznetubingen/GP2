#!/usr/bin/env python
import glob
import os

configfile: "config.yaml"

# Aligner
ALIGNER = config.get('aligner')

# INPUT FASTQ folder
FASTQDIR = config["resources"]["fastqdir"]

# BAM output
LRA_BAM = "{FASTQDIR}/{sample}_lra/alignment/{sample}_lra.bam"

# Input reference FASTA
FA_REF = config["resources"]["reference"]

# Reference index name
FA_REF_INDEX = FA_REF + ".gli"

# Parameter: sample_name
sample = config["resources"]["sample"]

# Parameter: threads

thread_n = config.get("thread_sv", 30)



######## RULES ##########        

# 1. Index REF   ##############################################################################

rule index_lra:
       input: 
           REF = FA_REF
       output:
           INDEX = FA_REF_INDEX
       conda: "../envs/SV_env.yml"
       threads: thread_n
       shell:
          "lra index -ONT {input.REF}"


# 2. Map with LRA and/or Winnowmap #############################################################  
  
rule map_lra:
       input:
          FQ = FASTQDIR,
          REF = FA_REF,
          INDEX = FA_REF_INDEX
       output:
          BAM = "{FASTQDIR}/{sample}_lra/alignment/{sample}_lra.bam",
          BAI = "{FASTQDIR}/{sample}_lra/alignment/{sample}_lra.bam.bai"
       conda: "../envs/SV_env.yml"
       threads: thread_n

       shell:
           "catfishq -r {input.FQ} | seqtk seq -A - | lra align -ONT -t {threads} {input.REF} - -p s | samtools addreplacerg -r \"@RG\tID:{sample}\tSM:{sample}\" - | samtools sort -@ {threads} -T {sample}_lra -O BAM -o {output.BAM} - && samtools index -@ {threads} {output.BAM}"


# 3. QC of alignment ############################################################################  

rule nanoplot_qc_lra:
        input:
            BAM = LRA_BAM
        output:
            DIR = directory("{FASTQDIR}/{sample}_lra/qc")
        params:
            sample = sample
        conda: "../envs/SV_env.yml"
        threads: thread_n
        shell:
            "NanoPlot -t {threads} --bam {input.BAM} --raw -o {output.DIR} -p {params.sample}_ --N50 --title {params.sample} --downsample 100000"



# 4. Call variants with CuteSV #######################################################################################

rule call_lra_cutesv:
        input:
            BAM = LRA_BAM,
            REF = FA_REF
        output:
            VCF_LRA = "{FASTQDIR}/{sample}_lra/sv_calls/{sample}_lra_cutesv_tmp.vcf"
        params:
            min_size = config.get("min_sv_length", 30),
            max_size = config.get("max_sv_length", 100000),
            min_read_support = 2,
            min_read_length = config.get("min_read_length", 1000),
            min_mq = config.get("min_read_mapping_quality", 20),
        conda: "../envs/SV_env.yml"
        threads: thread_n
        shell:
            "cuteSV -t {threads} --min_size {params.min_size} --max_size  {params.max_size} -S {sample} --retain_work_dir --report_readid --min_support {params.min_read_support} --genotype {input.BAM} {input.REF} {output.VCF_LRA} {sample}_lra/sv_calls/ "


# 5. Filter variants {scripts/wrapper_SV.py} ########################################################################

rule filter_lra_vcf:
        input:
             MOS = "{FASTQDIR}/{sample}_lra/depth",
             VCF = rules.call_lra_cutesv.VCF_LRA
        output:
             VCF = temp("{FASTQDIR}/{sample}_lra/sv_calls/{sample}_lra_cutesv_filtered_tmp.vcf")
        params:
            min_sv_length = config.get("min_sv_length", 30),
            max_sv_length = config.get("max_sv_length", 100000),
            target_bed = config.get("target_bed", None),
            sv_types = config.get("sv_type", "DEL INS INV DUP"),
        conda: "../envs/SV_env.yml"
        wrapper:
             f"file:scripts/wrapper_SV.py"

# 6. Sort and index VCF   ##############################################################################################

rule sort_lra_vcf:
        input:
             VCF = rules.filter_vcf_lra_cutesv.VCF
        output:
             VCF = temp("{FASTQDIR}/{sample}_lra/sv_calls/{sample}_cutesv_filtered.vcf")
        conda: "../envs/SV_env.yml"
        wrapper:
             "cat {input.VCF} | bgziptabix {output.VCF}"

    rule index_lra_vcf:
    input:
         VCF = rules.sort_lra_vcf.output.VCF
    output:
         VCF = "{FASTQDIR}/{sample}_lra/sv_calls/{sample}_cutesv_filtered.vcf.gz"
    conda: "../envs/SV_env.yml"
    shell:
         "cat {input.VCF} | bgziptabix {output.VCF}"









