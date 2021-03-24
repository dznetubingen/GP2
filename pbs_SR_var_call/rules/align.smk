#!/usr/bin/env python3

import glob
import os

configfile: "config.yaml"

FASTQDIR = config["resources"]["fastqdir"]
GENOMEDIR = config["resources"]["genomedir"]
SAMPLES, = glob_wildcards(FASTQDIR + "/{sample}_R1.fastq.gz")

#global singularity
singularity: "docker://continuumio/miniconda3"

checkpoint split_rg:
    input:
        R1_qc = "qc/{sample}_R1.fastq.gz",
        R2_qc = "qc/{sample}_R2.fastq.gz"
    output:
        rg_fq = temp(directory("rg_fq/{sample}"))
    container:
        "docker://quay.io/kmhernan/gdc-fastq-splitter"
    params:
        #gdc_fastq_splitter = config["tools"]["gdc_fastq_splitter"],
        prefix = "rg_fq/{sample}/{sample}_"
    threads: 1
#    resources:
#        mem = "1gb",
#        nodes = 1,
#        walltime = "12:00:00"
    shell:
        "mkdir -p rg_fq/{wildcards.sample};"
        "gdc-fastq-splitter -o {params.prefix} {input.R1_qc} {input.R2_qc}"

rule alignment:
    input:
        R1 = "rg_fq/{sample}/{sample}_{flowcell}_{lane}_R1.fq.gz",
        R2 = "rg_fq/{sample}/{sample}_{flowcell}_{lane}_R2.fq.gz",
        ref = config["resources"]["genome"]
    output:
        bam = temp("rg_fq/{sample}/aln/{sample}_{flowcell}_{lane}.bam")
    params:
        bwa_rg="@RG\\tID:{flowcell}.{lane}\\tCN:Macrogen\\tLB:{sample}\\tPL:illumina\\tPU:{flowcell}:{lane}\\tSM:{sample}"
    conda:
        "../envs/short_read_mapping.yml"
    threads: 15
#    resources:
#        mem= "1gb",
#        walltime = "24:00:00"
    shell:
        """
        bwa mem -M -R "{params.bwa_rg}" -t {threads} \
            {input.ref} {input.R1} {input.R2} \
            | samblaster -M \
            |samtools view -Sb - > {output.bam}
        """

rule sort_bam:
    input: rules.alignment.output
    output:
        bam = temp("rg_fq/{sample}/sort/{sample}_{flowcell}_{lane}.bam")
    conda:
        "../envs/short_read_mapping.yml"
    threads:
        10
#    resources:
#        mem = "5gb",
#        walltime = "04:00:00"
    shell:
        """
        samtools sort -t {threads} -T $TMPDIR -o {output.bam} {input}
        """

def alignment_input(wildcards):
    checkpoint_output = checkpoints.split_rg.get(sample = wildcards.sample).output[0]
    all_wildcards = glob_wildcards(os.path.join(checkpoint_output, "{sample}_{flowcell}_{lane}_R1.fq.gz"))
    all_files = []
    for sample, flowcell, lane in zip(all_wildcards.sample, all_wildcards.flowcell, all_wildcards.lane):
        all_files.append(f"rg_fq/{{sample}}/sort/{sample}_{flowcell}_{lane}.bam")
    return(all_files) 

rule merge_bam:
    input:
        alignment_input
    output:
        bam = temp("aln/{sample}_merge.bam")
    conda:
        "../envs/short_read_mapping.yml"
    threads:
        10
#    resources:
#        mem = "5gb",
#        walltime = "04:00:00"
    shell:   
        """
        mkdir -p aln
        samtools merge -@ {threads} - {input} | \
                samtools sort -@ {threads} -T $TMPDIR > {output.bam}
        """

rule set_tag:
    input:
        bam = "aln/{sample}_merge.bam",
        ref = config["resources"]["genome"]
    output:
        bam = temp("aln/{sample}_tag.bam")
    params:
        java_opts = "-Xmx8G -Dsamjdk.compression_level=5"
#    threads: 5
#    resources:
#        mem = "5gb",
#        walltime = "12:00:00"
    container:
        "docker://broadinstitute/gatk:4.1.9.0" 
    shell:
        """
         gatk --java-options "{params.java_opts}" \
               SetNmMdAndUqTags \
                      --INPUT {input.bam} \
                      --OUTPUT {output.bam} \
                      --REFERENCE_SEQUENCE {input.ref}
        """
# When using Markduplicates:
# Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
# This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
# While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
rule mark_duplicate:
    input: 
        "aln/{sample}_tag.bam"
    output:
        bam = "alignment/{sample}.bam", #protected("alignment/{sample}.bam")
        metrics = "stats/{sample}.metrics.txt"
    params:
        java_opts ="-Dsamjdk.compression_level=5 -Xms1g -Xmx48g -Djava.io.tmpdir=tmp_scratch",
        novaseq = "2500"
    threads: 1
#    envmodules:
#        "bio/gatk/4.1.2.0"
#    resources:
#        mem = "5gb",
#        walltime = "12:00:00" 
    container:
        "docker://broadinstitute/gatk:4.1.9.0" 
    shell:
        """
        gatk --java-options "{params.java_opts}" \
             MarkDuplicates \
             --INPUT {input} \
             --OUTPUT {output.bam} \
             --METRICS_FILE {output.metrics} \
             --VALIDATION_STRINGENCY SILENT \
             --OPTICAL_DUPLICATE_PIXEL_DISTANCE {params.novaseq} \
             --ASSUME_SORT_ORDER queryname \
             --CREATE_MD5_FILE true \
             --CREATE_INDEX true
        """
###For running MarkDuplicatesSpark 
#rule mark_duplicate:
#    input:
#        "aln/{sample}_tag.bam"
#    output:
#        bam = protected("alignment/{sample}.bam"),
#        metrics = "stats/{sample}.metrics.txt"
#    params:
#        java_opts ="-Dsamjdk.compression_level=5 -Xms48G -Djava.io.tmpdir=tmp_scratch",
#        novaseq_dup= "2500"
#    threads: 5
#    resources:
#        mem = "5gb",
#        walltime = "12:00:00"
#    container:
#        "docker://broadinstitute/gatk:4.1.9.0"
#    shell:
#        """
#        gatk --java-options "{params.java_opts}" \
#            MarkDuplicatesSpark  \
#                -I {input} \
#                -O {output.bam} \
#                --optical-duplicate-pixel-distance {params.novaseq_dup} \
#                -OBI true \
#                -ASSUME_SORT_ORDER queryname \
#                --CREATE_MD5_FILE true \
#                --CREATE_INDEX true \
#                --conf 'spark.executor.cores={threads}' \
#                --conf 'spark.local.dir={params.tmpdir}'
#        """


rule aln_stat:
    input:
        "alignment/{sample}.bam"
    output:
        "stats/{sample}_aln.stats"
    conda:
        "../envs/short_read_mapping.yml"
    #params:
    #    sambamba= config["tools"]["sambamba"]
    threads:
        5
    shell:
        """
        sambamba flagstat -t {threads} {input} > {output}
    
        """

