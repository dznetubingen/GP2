# This snakemake pipeline is for the annotation and filter for candidate variants

**Input: The output (vcf file) from jointgenotyping workflow (Terra)**

## Steps:
**Prepare variants for filtering to obtain candidate causal vairants:**
[commands](https://github.com/dznetubingen/GP2/blob/d71ac897f1ce9b195279fb06e43c4e4388df4eed/Snakefile)
1. Normalize indels and split multialleic variants (see [explanation](https://genome.sph.umich.edu/wiki/Variant_Normalization))
2. Split the vcf by chromosome for parallel processing 
3. For the vairants on each chromosome: 
- annotate the vcf with vep and split the vep annotation with [bcftools +split-vep](https://samtools.github.io/bcftools/howtos/plugin.split-vep.html) 
- calculate CADD scores
- calculate genotype counts in GP2 cohort 
4. Add CADD scores, genotype counts for each variant in GP2 cohort and AMP-PD cases, topmed freq as well as gnomAD freq and nhomalt (number of alternative homozygous) to vcf INFO for downstream filtering

**Filtering** (here [slivar](https://github.com/brentp/slivar) is used):
Global filter by gnomAD_nhomalt < 10 and gnomAD_af < 0.05

For each family (see [rule slivar_filter](https://github.com/dznetubingen/GP2/blob/710b391b1b57d3f335441d46370b5c1b3f8df7b2/rules/slivar_filter.smk) for the commands):
- (segragating) recessive (gnomAD_af < 0.05, the affected individuals are alternative homozygous)
- (segragating) dominant (gnomAD_af < 0.01, the affected individuals are heterozygous)
- x-reccessive (same as segratating recessive)
- x-denovo  (x-dominant)
- compound-heterozygotes (gnomAD_nhomalt < 10)

** Check [slivar-function.js](https://github.com/dznetubingen/GP2/blob/710b391b1b57d3f335441d46370b5c1b3f8df7b2/slivar-functions.js) for the function expression


[Workflow visulization](https://github.com/dznetubingen/GP2/blob/710b391b1b57d3f335441d46370b5c1b3f8df7b2/rulegraph.pdf)


