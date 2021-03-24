#!/usr/bin/env python3

import glob
import os

configfile: "config.yaml"
FASTQDIR = config["resources"]["fastqdir"]
GENOMEDIR = config["resources"]["genomedir"]

##Parsing wildcards
chr_list = list(range(1,22))+["X", "Y", "M"]
chr_id= ["chr" + str(i) for i in chr_list]
#int_list=" ".join (   [f"-L {mychr}" for mychr in chr_list   ]     )
#print (f"Number of samples: {len (sample)}" )
#print (sample)

SAMPLES, = glob_wildcards(FASTQDIR + "/{sample}_R1.fastq.gz")


rule gatk_baserecalibrator:
    input:
        bam = "alignment/{sample}.bam",
        ref = config["resources"]["genome"] ,
        dbsnp =  GENOMEDIR + "/Homo_sapiens_assembly38.dbsnp138.vcf",
        known_indel = GENOMEDIR + "/Homo_sapiens_assembly38.known_indels.vcf.gz",
        interval = GENOMEDIR + "/wgs_calling_regions.hg38.interval_list"
    output:
        "bam_recal/{sample}_recal.table"
    params:
        java_opts = "-Xmx8G -Djava.io.tmpdir=tmp_scratch"
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
#    resources:
#        mem = "8gb",
#        walltime = "12:00:00"
#    benchmark:
#        "benchmarks/recal_{sample}.tsv"
    shell:
        """
         gatk --java-options "{params.java_opts}" \
              BaseRecalibrator \
              -R {input.ref} \
              -I {input.bam} \
              --use-original-qualities \
              -O {output} \
              --known-sites {input.dbsnp} \
              --known-sites {input.known_indel} \
              -L {input.interval}
        """

rule apply_BQSR:
    input:
        recal = rules.gatk_baserecalibrator.output,
        bam = "alignment/{sample}.bam",
        ref = config["resources"]["genome"],
        interval = GENOMEDIR + "/wgs_calling_regions.hg38.interval_list"
    output:
        "bam_recal/{sample}_recal.bam"
    params:
        java_opts = "-Xmx12G -Dsamjdk.compression_level=5 -Djava.io.tmpdir=tmp_scratch"
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
#    resources:
#        mem = "12gb",
#        walltime = "24:00:00"
#    benchmark:
#        "benchmarks/bqsr_{sample}.tsv"
    shell:
        """
         gatk --java-options "{params.java_opts}" \
            ApplyBQSR \
            -R {input.ref} \
            -I {input.bam} \
            -O {output} \
            -L {input.interval} \
            --bqsr-recal-file {input.recal} \
            --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
            --add-output-sam-program-record \
            --create-output-bam-md5 true \
            --use-original-qualities
        """

rule index:
    input:
        bam = "bam_recal/{sample}_recal.bam"
    output:
        "bam_recal/{sample}_recal.bam.bai"
    conda:
        "../envs/short_read_mapping.yml"
    shell:
        """
        samtools index {input.bam}
        """

rule haplotype_caller:
    input:
        recal_bam = rules.apply_BQSR.output,
        ref = config["resources"]["genome"]
    output:
        "gvcf/{chromosome}/{sample}_{chromosome}.g.vcf.gz"
    params:
        chr = "{chromosome}",
        java_opts = "-Xmx12G -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=tmp_scratch"    
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
    shell:
        """
        gatk --java-options "{params.java_opts}" HaplotypeCaller \
         -I {input.recal_bam} \
         -R {input.ref} \
         -L {params.chr} \
         -O {output} \
         -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
         -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
         --ERC GVCF
        """

# Gather all samples given chromosome

rule GenomicsDB_import:
    input:
        expand("gvcf/{{chromosome}}/{sample}_{{chromosome}}.g.vcf.gz", sample=SAMPLES)
    output:
        samp_map = "db_imp/{chromosome}.map",
        db_dir = temp(directory("db_imp/db_{chromosome}"))
    params:
        chr = "{chromosome}",
        ref = config["resources"]["genome"]
        java_opts = "-Xmx4g -Xms4g"
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
    shell:
        " paste -- " +
        " <(echo {input}|xargs -n 1 basename -s \"_{params.chr}.g.vcf.gz\") " +
        " <(echo {input}|tr ' ' '\\n' ) > {output.samp_map} \n " +
        """
        gatk --java-options "{params.java_opts}" \
             GenomicsDBImport \
             --sample-name-map {output.samp_map} \
             --genomicsdb-workspace-path {output.db_dir} \
             -L {params.chr} \
             -R {params.ref} \
             --tmp-dir=./tmp_scratch 
             --reader-threads 5
        """

rule Joint_Genotyping:
    input:
        "db_imp/db_{chromosome}"
    output:
        "raw_vcf/{chromosome}.vcf.gz"
    params:
        chr = "{chromosome}",
        ref = config["resources"]["genome"],
        dbsnp =  GENOMEDIR + "/Homo_sapiens_assembly38.dbsnp138.vcf",
        java_opts = "-Xmx5g -Xms5g"
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
    shell:
        """
        gatk --java-options "{params.java_opts}" \
             GenotypeGVCFs \
             -R {params.ref} \
             -L {params.chr} \
             -D {params.dbsnp} \
             -O {output} \
             -G StandardAnnotation \
             --only-output-calls-starting-in-intervals \
             --use-new-qual-calculator \
             -V gendb://{input}
        """
