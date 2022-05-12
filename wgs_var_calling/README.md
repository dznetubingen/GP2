# The Short-read variant calling pipelines for GP2 samples on Terra

**Input: Paired-end fastq files**

- 1-SplitRG_Paired_Fastq.wdl: Split paried-end fastq files by read groups.
- 2-Paired-fastq-to-unmapped-bam.wdl: Convert fastqs into ubam files
- 3-PairedEndSingleSampleWf.wdl: Generate individual cram, gvcf and WGS metrics from ubams
- 4-generate-sample-map.wdl: Generate a list of gvcfs used for joint genotyping
- 5-joint-discovery-gatk4-gp2.wdl: Perfrom joint genotyping across the samples to produce a callset. 
