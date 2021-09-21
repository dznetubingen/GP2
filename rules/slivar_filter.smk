#!/usr/bin/env python
import glob
from pathlib import Path
from snakemake.utils import min_version

shell.prefix("set -o pipefail; umask 002; ")  # set g+w

if not os.path.exists("short_list"):
    os.makedirs("short_list")

rule slivar_filter:
    input:
        vcf = "final_vcfs/{chromosome}.vcf.gz",
        ped = "slivar/GP2_trimmed.ped"
    output:
        "slivar_filtered/{chromosome}.vcf"
    params:
        chr = "{chromosome}"
    container:
        "docker://brentp/slivar:v0.2.5"
    resources:
        nodes = 1
    shell:
        """
        base=$(basename {input.vcf})
        if [ $base == "chrX.vcf.gz" ];
        then
           slivar expr --vcf {input.vcf} \
                       --ped {input.ped} \
                       -o {output} \
                       --js slivar/slivar-functions.js \
                       --pass-only \
                       --info 'INFO.gnomad_nhomalt < 10 && variant.ALT[0] != "*" && variant.call_rate > 0.95' \
                       --family-expr 'x_denovo:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af < 0.01' \
                       --family-expr 'x_recessive:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_recessive_x) && INFO.gnomad_popmax_af < 0.05' \
                       --trio 'comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_nhomalt < 10'
          else 
             slivar expr --vcf {input.vcf} \
                         --ped {input.ped} \
                         -o {output} \
                         --js slivar/slivar-functions.js \
                         --pass-only \
                         --info 'INFO.gnomad_nhomalt < 10 && variant.ALT[0] != "*" && variant.call_rate > 0.95' \
                         --family-expr 'recessive:fam.every(segregating_recessive) && INFO.gnomad_popmax_af < 0.05' \
                         --family-expr 'x_denovo:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af < 0.01' \
                         --family-expr 'x_recessive:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_recessive_x) && INFO.gnomad_popmax_af < 0.05' \
                         --family-expr 'dominant:INFO.gnomad_popmax_af < 0.01 && fam.every(segregating_dominant)' \
                         --trio 'comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_nhomalt < 10' \
                         --family-expr 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.01'
          fi
          """

rule slivar_comphet:
    input:
        vcf = "slivar_filtered/{chromosome}.vcf",
        ped = "slivar/GP2_trimmed.ped"
    output:
        "slivar_filtered/comphet/{chromosome}.vcf"
    params:
        chr = "{chromosome}"
    container:
        "docker://brentp/slivar:v0.2.5"
    resources:
        nodes = 1
    shell:
        """
        slivar compound-hets -v {input.vcf} \
               --allow-non-trios --sample-field comphet_side --sample-field denovo -p {input.ped} \
               -o {output}
        """

rule slivar_tsv:
    input:
        vcf = rules.slivar_filter.output,
        ped = "slivar/GP2_trimmed.ped"
    output:
        "short_list/{chromosome}.tsv"
    params:
        chr = "{chromosome}"
    container:
        "docker://brentp/slivar:v0.2.5"
    resources:
        nodes = 1
    shell:
        """
        base=$(basename {input.vcf})
        if [ $base == "chrX.vcf" ];
        then
           slivar tsv \
                  -s x_denovo \
                  -s x_recessive \
                  -i gnomad_popmax_af -i gnomad_nhomalt -i topmed_af -i CADD_RAW -i CADD_PHRED -i GP2_nhet -i GP2_nhomalt -i GP2_nmiss -i impactful -i genic \
                  -c CSQ \
                  --csq-column Existing_variation \
                  --csq-column gnomAD_genomes_controls_and_biobanks_nhomalt \
                  --csq-column BIOTYPE \
                  --csq-column MAX_AF_POPS \
                  --csq-column clinvar_clnsig \
                  --csq-column clinvar_trait \
                  -p {input.ped} \
                  {input.vcf} > {output}
        else 
            slivar tsv \
                   -s x_denovo \
                   -s recessive \
                   -s x_recessive \
                   -s dominant \
                   -i gnomad_popmax_af -i gnomad_nhomalt -i topmed_af -i CADD_RAW -i CADD_PHRED -i GP2_nhet -i GP2_nhomalt -i GP2_nmiss -i impactful -i genic \
                   -c CSQ \
                   --csq-column Existing_variation \
                   --csq-column gnomAD_genomes_controls_and_biobanks_nhomalt \
                   --csq-column BIOTYPE \
                   --csq-column MAX_AF_POPS \
                   --csq-column clinvar_clnsig \
                   --csq-column clinvar_trait \
                   -p {input.ped} \
                  {input.vcf} > {output}
        fi 
        """

rule comphet_tsv:
    input:
        vcf = rules.slivar_comphet.output,
        ped = "slivar/GP2_trimmed.ped"
    output:
        "short_list/comphet_{chromosome}.tsv"
    params:
        chr = "{chromosome}"
    container:
        "docker://brentp/slivar:v0.2.5"
    resources:
        nodes = 1
    shell:
        """
        slivar tsv \
            -s slivar_comphet \
            -i gnomad_popmax_af -i gnomad_nhomalt -i topmed_af -i CADD_RAW -i CADD_PHRED -i GP2_nhet -i GP2_nhomalt -i GP2_nmiss -i impactful -i genic \
             -c CSQ \
             --csq-column Existing_variation \
             --csq-column gnomAD_genomes_controls_and_biobanks_nhomalt \
             --csq-column BIOTYPE \
             --csq-column MAX_AF_POPS \
             --csq-column clinvar_clnsig \
             --csq-column clinvar_trait \
             -p {input.ped} \
             {input.vcf} \
             | {{ grep -v ^# || true; }} >> {output}
        """


rule merge_tsv:
    input: 
        mode= lambda wildcards: expand("short_list/{chromosome}.tsv", chromosome= chr_id),
        comphet = lambda wildcards: expand("short_list/comphet_{chromosome}.tsv", chromosome= chr_id)
    output:
        "short_list/all_chrs.tsv"
    resources:
        nodes = 1
    shell:
        """
        awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input.mode} {input.comphet} > {output}
        """
