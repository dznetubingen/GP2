#!/bin/env python
import glob
import os

configfile: "config.yaml"
FASTQDIR = config["resources"]["fastqdir"]
GENOMEDIR = config["resources"]["genomedir"]

##Parsing wildcards
chr_list = list(range(1,22))+["X", "Y","M"]
chr_id= ["chr" + str(i) for i in chr_list]
#int_list=" ".join (   [f"-L {mychr}" for mychr in chr_list   ]     )
#print (f"Number of samples: {len (sample)}" )
#print (sample)

SAMPLES, = glob_wildcards(FASTQDIR + "/{sample}_R1.fastq.gz")

rule gather_rawvcf:
    input:
        lambda w: expand("raw_vcf/{chromosome}.vcf.gz", chromosome= chr_id)
        #expand("raw_vcf/{chromosome}.vcf.gz", chromosome= chr_id)
    output:
        temp("VQSR_all/FD_10L_raw.vcf.gz")
    params:
        java_opts = "-Xmx3g -Xms3g",
        vcf = lambda w: "-I=" + " -I=".join(expand("raw_vcf/{chromosome}.vcf.gz", chromosome= chr_id))
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
#    resources:
#        mem = "12gb",
#        walltime = "24:00:00"
#    benchmark:
#        "benchmarks/bqsr_{sample}.tsv"
    shell:
        """
        gatk  --java-options "{params.java_opts}" \
              GatherVcfs \
              {params.vcf} \
              -O {output} && tabix -p vcf {output}
        """

rule make_site:
    input:
       "VQSR_all/FD_10L_raw.vcf.gz"        
    output:
        temp("VQSR_all/FD_10L_sitesonly.vcf.gz")
    params:
        java_opts = "-Xmx3g -Xms3g"
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
#    resources:
#        mem = "12gb",
#        walltime = "24:00:00"
#    benchmark:
#        "benchmarks/bqsr_{sample}.tsv"
    shell:
        """
        gatk  --java-options "{params.java_opts}" \
              MakeSitesOnlyVcf \
              -I {input} \
              -O {output}
        """

rule VQSR_indel:
    input:
        rules.make_site.output
    output:
        recal = temp("VQSR_all/indel/FD_10L.recal"),
        tranches = temp("VQSR_all/indel/FD_10L.tranches")
    params:
        mills_vcf = GENOMEDIR + "/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        axiomPoly_vcf = GENOMEDIR + "/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz",
        dbsnp_vcf = GENOMEDIR + "/Homo_sapiens_assembly38.dbsnp138.vcf",
        java_opts = "-Xmx24g -Xms24g"    
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
    shell:
        """
        gatk --java-options "{params.java_opts}" \
               VariantRecalibrator \
               -V {input}  \
               -O {output.recal}  \
               --tranches-file {output.tranches} \
               --trust-all-polymorphic \
               -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
               -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
               -mode INDEL \
               --max-gaussians 4 \
               --resource:mills,known=false,training=true,truth=true,prior=12 {params.mills_vcf} \
               --resource:axiomPoly,known=false,training=true,truth=false,prior=10 {params.axiomPoly_vcf} \
               --resource:dbsnp,known=true,training=false,truth=false,prior=2 {params.dbsnp_vcf}
        """

rule VQSR_snp:
    input:
        rules.make_site.output
    output:
        recal = temp("VQSR_all/snp/FD_10L.recal"),
        tranches = temp("VQSR_all/snp/FD_10L.tranches")
    params:
        hapmap_vcf = GENOMEDIR + "/hapmap_3.3.hg38.vcf.gz",
        omni_vcf = GENOMEDIR + "/1000G_omni2.5.hg38.vcf.gz",
        thousandG_vcf = GENOMEDIR + "/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp_vcf = GENOMEDIR + "/Homo_sapiens_assembly38.dbsnp138.vcf",
        java_opts = "-Xmx50g -Xms50g"
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
    shell:
        """
        gatk --java-options "{params.java_opts}" \
              VariantRecalibrator \
              -V {input} \
              -O {output.recal} \
              --tranches-file {output.tranches} \
              --trust-all-polymorphic \
              -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
              -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
              -mode SNP \
              --max-gaussians 6 \
              --resource:hapmap,known=false,training=true,truth=true,prior=15 {params.hapmap_vcf} \
              --resource:omni,known=false,training=true,truth=true,prior=12 {params.omni_vcf} \
              --resource:1000G,known=false,training=true,truth=false,prior=10 {params.thousandG_vcf} \
              --resource:dbsnp,known=true,training=false,truth=false,prior=7 {params.dbsnp_vcf}
        """

rule filter_indel:
    input:
        vcf = "VQSR_all/FD_10L_raw.vcf.gz",
        recal = "VQSR_all/indel/FD_10L.recal",
        tranches = "VQSR_all/indel/FD_10L.tranches"
    output:
        temp("filtered_vcf/FD_10L_indel.vcf.gz")
    params:
        java_opts = "-Xmx5g -Xms5g"
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
    shell:
        """
        gatk --java-options "{params.java_opts}" \
             ApplyVQSR \
             -V {input.vcf} \
             --recal-file {input.recal} \
             --tranches-file {input.tranches} \
             --truth-sensitivity-filter-level 99.7 \
             --create-output-variant-index true \
             -mode INDEL \
             -O {output}
        """

rule filter_snp:
    input:
        vcf = "filtered_vcf/FD_10L_indel.vcf.gz",
        recal = "VQSR_all/snp/FD_10L.recal",
        tranches = "VQSR_all/snp/FD_10L.tranches"
    output:
        "filtered_vcf/FD_10L.vcf.gz"
    params:
        java_opts = "-Xmx5g -Xms5g"
    container:
        "docker://broadinstitute/gatk:4.1.9.0"
    shell:
        """
        gatk --java-options "{params.java_opts}" \
             ApplyVQSR \
             -V {input.vcf} \
             --recal-file {input.recal} \
             --tranches-file {input.tranches} \
             --truth-sensitivity-filter-level 99.7 \
             --create-output-variant-index true \
             -mode SNP \
             -O {output} \
        """

