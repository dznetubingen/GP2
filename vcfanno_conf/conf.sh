#!/bin/bash
for i in {1..22} X Y

do
echo -e "[[annotation]]\nfile=\"CADD/chr${i}.tsv.gz\"\nnames=[\"CADD_RAW\", \"CADD_PHRED\"]\nops=[\"self\", \"self\"]\ncolumns=[5, 6]\n
[[annotation]]\nfile=\"vcfs/gcount/chr${i}.vcf.gz\"\nfields = [\"Cases\", \"Controls\"]\nnames=[\"GP2_affected_gt\", \"GP2_unaffected_gt\"]\nops=[\"self\", \"self\"]\n
[[annotation]]\nfile=\"gcounts/amp-pd_gcount/chr${i}.gcount.gz\"\nnames=[\"amp-pd_case_nhet-total\", \"amp-pd_case_nhomalt-total\",\"amp-pd_case_nmiss-total\"]\nops=[\"self\", \"self\",\"self\"]\ncolumns=[7,8,9]\n
[[annotation]]\nfile=\"/mnt/vol009/hg38_vep105/slivar/CNCR_scores_hits.bed.gz\"\ncolumns=[4]\nnames=[\"CNCR\"]\nops=[\"mean\"]" > vcfanno_chr${i}.conf
done
