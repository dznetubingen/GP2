#!/bin/bash
for i in {1..22} X Y

do
echo -e "[[annotation]]\nfile=\"CADD/chr${i}.tsv.gz\"\nnames=[\"CADD_RAW\", \"CADD_PHRED\"]\nops=[\"self\", \"self\"]\ncolumns=[5, 6]\n
[[annotation]]\nfile=\"gp2_gcount/chr${i}.gcount.gz\"\nnames=[\"GP2_nhet\", \"GP2_nhomalt\",\"GP2_nmiss\"]\nops=[\"self\", \"self\",\"self\"]\ncolumns=[7,8,9]\n
[[annotation]]\nfile=\"amp-pd_gcount/chr${i}.gcount.gz\"\nnames=[\"amp-pd_case_nhet\", \"amp-pd_case_nhomalt\",\"amp-pd_case_nmiss\"]\nops=[\"self\", \"self\",\"self\"]\ncolumns=[7,8,9]" > vcfanno_chr${i}.conf
done
