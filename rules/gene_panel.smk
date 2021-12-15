#!/usr/bin/env python
import glob
from pathlib import Path
from snakemake.utils import min_version

shell.prefix("set -o pipefail; umask 002; ")  # set g+w

if not os.path.exists("short_list"):
    os.makedirs("short_list")

rule filter_frq:
    input:
        vcf = "final_vcfs/GP2_annotate.vcf.gz",
        bed = "PD_panel.bed"   #include both core and extensive panels
    output:
        "panel/GP2_annotate_filtered.vcf.gz"
    container:
        "docker://brentp/slivar:v0.2.7"
    resources:
        nodes = 1
    shell:
        """
        slivar expr --vcf {input} \
                    -o {output} \
                    --js slivar/slivar-functions.js \
                    --pass-only \
                    --info 'INFO.gnomad_popmax_af < 0.05 && variant.ALT[0] != "*" && variant.call_rate > 0.95' \
                    --region {input.bed}
        """

rule vcf2db:
    input:
        rules.filter_frq.output
    output:
        "panel/gp2_panel.db"
    conda:
        "../envs/vcf2db.yml"
    resources:
        nodes = 1
    shell:
        """
        vcf2db.py --expand gt_quals --expand gt_depths --expand gt_alt_depths --expand gt_ref_depths --expand gt_types \
                  --a-ok MLEAC --a-ok MLEAF --a-ok AF --a-ok AC --a-ok AN \
                  --a-ok gnomad_af --a-ok gnomad_ac --a-ok gnomad_nhomalt --a-ok topmed_af --a-ok GP2_nhet --a-ok GP2_nhomalt --a-ok GP2_nmiss --a-ok CADD_RAW --a-ok CADD_PHRED --a-ok amp-pd_case_nhet --a-ok amp-pd_case_nhomalt --a-ok amp-pd_case_nmiss \
{input.vcf} {config[gp2_ped]} {output}
        """
        
         
  

rule gemini_query:



