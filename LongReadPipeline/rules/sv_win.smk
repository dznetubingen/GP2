#!/usr/bin/env python
import glob
import os

configfile: "config.yaml"


# INPUT FASTQ folder
FASTQDIR = config["resources"]["fastqdir"]

# BAM output
WIN_BAM = "{FASTQDIR}/{sample}_win/alignment/{sample}_winnowmap.bam"



# Input reference FASTA
FA_REF = config["resources"]["reference"]

# Parameter: sample_name
sample = config["resources"]["sample"]

# Parameter: threads

thread_n = config.get("thread_sv", 30)



######## RULES ##########        

# 1. Index REF   ##############################################################################

rule kmer_winnowmap:
       input:
           INPUT = FASTQDIR         
           REF = FA_REF
       output:
           REP15 = "{FASTQDIR}/{sample}_repetitive_k15.txt"
       threads: thread_n
       shell:
          """
          meryl count k=15 output {input.INPUT}/{sample}_merylDB {input.REF}
          meryl print greater-than distinct=0.9998 {input.INPUT}/{sample}_merylDB > {output}
          """


# 2. Map with LRA and/or Winnowmap #############################################################  
  
rule map_winnowmap:
       input:
          FQ = FASTQDIR,
          REF = FA_REF,
          REP15 = rules.kmer_winnowmap.output.REP15,
       output:
          BAM = "{FASTQDIR}/{sample}_win/alignment/{sample}_winnowmap.bam",
          BAI = "{FASTQDIR}/{sample}_win/alignment/{sample}_winnowmap.bam.bai"
       singularity: "docker://avalillarionova/lrp_tools:latest"
       threads: thread_n
       shell:
           """
             conda activate SV-env
             catfishq -r {input.FQ} | seqtk seq -A - | winnowmap -W {input.REP15} -ax map-ont {input.REF} - | samtools addreplacerg -r \"@RG\tID:{sample}\tSM:{sample}\" - | samtools sort -@ {threads} -T {sample}_win -O BAM -o {output.BAM} - && samtools index -@ {threads} {output.BAM}
           """


# 3. QC of alignment ############################################################################  

rule nanoplot_qc_winnowmap:
        input:
            BAM = WIN_BAM
        output:
            DIR = directory("{FASTQDIR}/{sample}_win/qc")
        params:
            sample = sample
        conda: "../envs/SV_env.yml"
        threads: thread_n
        shell:
            "NanoPlot -t {threads} --bam {input.BAM} --raw -o {output.DIR} -p {params.sample}_ --N50 --title {params.sample} --downsample 100000"



# 4. Call variants with CuteSV #######################################################################################

rule call_win_cutesv:
        input:
            BAM = WIN_BAM,
            REF = FA_REF
        output:
            VCF_WIN = "{FASTQDIR}/{sample}_win/sv_calls/{sample}_win_cutesv_tmp.vcf"
        params:
            min_size = config.get("min_sv_length", 30),
            max_size = config.get("max_sv_length", 100000),
            min_read_support = 2,
            min_read_length = config.get("min_read_length", 1000),
            min_mq = config.get("min_read_mapping_quality", 20),
        conda: "../envs/SV_env.yml"
        threads: thread_n
        shell:
            "cuteSV -t {threads} --min_size {params.min_size} --max_size  {params.max_size} -S {sample} --retain_work_dir --report_readid --min_support {params.min_read_support} --genotype {input.BAM} {input.REF} {output.VCF_WIN} {sample}_win/sv_calls/ "


# 5. Filter variants {scripts/wrapper_SV.py} ########################################################################

rule filter_win_vcf:
        input:
             MOS = "{FASTQDIR}/{sample}_win/depth",
             VCF = rules.call_win_cutesv.VCF_WIN
        output:
             VCF = temp("{FASTQDIR}/{sample}_win/sv_calls/{sample}_win_cutesv_filtered_tmp.vcf")
        params:
            min_sv_length = config.get("min_sv_length", 30),
            max_sv_length = config.get("max_sv_length", 100000),
            target_bed = config.get("target_bed", None),
            sv_types = config.get("sv_type", "DEL INS INV DUP"),
        conda: "../envs/SV_env.yml"
        wrapper:
             f"file:scripts/wrapper_SV.py"

# 6. Sort and index VCF   ##############################################################################################

rule sort_win_vcf:
        input:
             VCF = rules.filter_win_lra_cutesv.VCF
        output:
             VCF = temp("{FASTQDIR}/{sample}_win/sv_calls/{sample}_cutesv_filtered.vcf")
        conda: "../envs/SV_env.yml"
        wrapper:
             "cat {input.VCF} | bgziptabix {output.VCF}"

    rule index_lra_vcf:
    input:
         VCF = rules.sort_win_vcf.output.VCF
    output:
         VCF = "{FASTQDIR}/{sample}_win/sv_calls/{sample}_cutesv_filtered.vcf.gz"
    conda: "../envs/SV_env.yml"
    shell:
         "cat {input.VCF} | bgziptabix {output.VCF}"









