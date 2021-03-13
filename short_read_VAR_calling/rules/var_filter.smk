#!/bin/env python
import glob
import os

if not os.path.exists("tmp_scratch"):
    os.makedirs("tmp_scratch")

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

rule ExcessHet:
    input:
        vcf = "raw_vcf/{chromosome}.vcf.gz"
    output:
        vcf="cohort_excesshet_{chromosome}.vcf.gz" if config["condition"] else []
    params:
        java_opts = "-Xmx3g -Xms3g"
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
#    resources:
#        mem = "8gb",
#        walltime = "12:00:00"
#    benchmark:
#        "benchmarks/recal_{sample}.tsv"
    shell:
        """
         num=$(ls {input} | wc -l)
         if [ "$num" -gt 100]
              then
              gatk --java-options "{params.java_opts}" VariantFiltration \
                    -V cohort.vcf.gz \
                    --filter-expression "ExcessHet > 54.69" \
                    --filter-name ExcessHet \
                    -O cohort_excesshet.vcf.gz 
         fi
        """

rule make_site:
    input:
        
        ref = ["resources"]["genome"],
    output:
        bam = "bam_recal/{sample}_recal.bam"
    params:
        java_opts = "-Xmx12G -Dsamjdk.compression_level=5"
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
#    resources:
#        mem = "12gb",
#        walltime = "24:00:00"
#    benchmark:
#        "benchmarks/bqsr_{sample}.tsv"
    shell:
        """
         gatk MakeSitesOnlyVcf \
               -I cohort_excesshet.vcf.gz \
               -O cohort_sitesonly.vcf.gz
        """
rule VQSRD_indel:
    input:
        recal_bam = rules.apply_BQSR.output,
        ref = ["resources"]["genome"]
    output:
        "gvcf/{chromosome}/{sample}_{chromosome}.g.vcf.gz"
    params:
        chr = "{chromosome}",
        java_opts = "-Xmx12G -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"    
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
    shell:
        """
        gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
               -V cohort_sitesonly.vcf.gz \
               --trust-all-polymorphic \
               -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
               -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \      
               -mode INDEL \
               --max-gaussians 4 \
               -resource:mills,known=false,training=true,truth=true,prior=12:Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
               -resource:axiomPoly,known=false,training=true,truth=false,prior=10:Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
               -resource:dbsnp,known=true,training=false,truth=false,prior=2:Homo_sapiens_assembly38.dbsnp138.vcf \
               -O cohort_indels.recal \
               --tranches-file cohort_indels.tranches
        """

# Gather all samples given chromosome

rule VQSRD_snp:
    input:
        expand("gvcf/{{chromosome}}/{sample}_{{chromosome}}.g.vcf.gz", sample=SAMPLES)
    output:
        samp_map = "db_imp/{chromosome}.map",
        db_dir = temp(directory("db_imp/db_{chromosome}"))
    params:
        chr = "{chromosome}",
        ref = ["resources"]["genome"],
        java_opts = "-Xmx4g -Xms4g"
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
    shell:
        """
        gatk --java-options "-Xmx3g -Xms3g" VariantRecalibrator \
              -V cohort_sitesonly.vcf.gz \
              --trust-all-polymorphic \
              -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
              -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
              -mode SNP \
              --max-gaussians 6 \
              -resource:hapmap,known=false,training=true,truth=true,prior=15:hapmap_3.3.hg38.vcf.gz \
              -resource:omni,known=false,training=true,truth=true,prior=12:1000G_omni2.5.hg38.vcf.gz \
              -resource:1000G,known=false,training=true,truth=false,prior=10:1000G_phase1.snps.high_confidence.hg38.vcf.gz \
              -resource:dbsnp,known=true,training=false,truth=false,prior=7:Homo_sapiens_assembly38.dbsnp138.vcf \
              -O cohort_snps.recal \
              --tranches-file cohort_snps.tranches
        """

rule filter_indel:
    input:
        "db_imp/db_{chromosome}"
    output:
        "vcf/{chromosome}_gatk4.vcf.gz"
    params:
        chr = "{chromosome}",
        ref = ["resources"]["genome"]
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
    shell:
        """
        gatk --java-options "-Xmx5g -Xms5g" \
             ApplyVQSR \
             -V cohort_excesshet.vcf.gz \
             --recal-file cohort_indels.recal \
             --tranches-file cohort_indels.tranches \
             --truth-sensitivity-filter-level 99.7 \
             --create-output-variant-index true \
             -mode INDEL \
             -O indel.recalibrated.vcf.gz
        """
rule filter_snp:
    input:
        "db_imp/db_{chromosome}"
    output:
        "vcf/{chromosome}_gatk4.vcf.gz"
    params:
        chr = "{chromosome}",
        ref = ["resources"]["genome"]
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
    shell:
        """
        gatk --java-options "-Xmx5g -Xms5g" \
             ApplyVQSR \
             -V indel.recalibrated.vcf.gz \
             --recal-file ${snps_recalibration} \
             --tranches-file ${snps_tranches} \
             --truth-sensitivity-filter-level 99.7 \
             --create-output-variant-index true \
             -mode SNP \
             -O snp.recalibrated.vcf.gz 
        """

