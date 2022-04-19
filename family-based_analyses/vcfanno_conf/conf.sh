#!/bin/bash
for i in {1..22} X Y

do
echo -e "[[annotation]]\nfile=\"CADD/chr${i}.tsv.gz\"\nnames=[\"CADD_RAW\", \"CADD_PHRED\"]\nops=[\"self\", \"self\"]\ncolumns=[5, 6]\n
[[annotation]]\nfile=\"vcfs/gcount/chr${i}.vcf.gz\"\nfields = [\"Cases\", \"Controls\"]\nnames=[\"GP2_affected_gt\", \"GP2_unaffected_gt\"]\nops=[\"self\", \"self\"]\n
[[annotation]]\nfile=\"amp-pd-gcounts/chr${i}.vcf.gz\"\nfields = [\"Cases\"]\nnames=[\"amp_pd_cases_gt\"]\nops=[\"self\"]\n
[[annotation]]\nfile=\"/mnt/vol009/hg38_vep105/slivar/CNCR_scores_hits.bed.gz\"\ncolumns=[4]\nnames=[\"CNCR\"]\nops=[\"mean\"]" > vcfanno_chr${i}.conf
done



#"\nnames=[\"amp-pd_case_nhet_3359\", \"amp-pd_case_nhomalt_3359\",\"amp-pd_case_nmiss_3359\"]\nops=[\"self\", \"self\",\"self\"]\ncolumns=[7,8,9]\n
