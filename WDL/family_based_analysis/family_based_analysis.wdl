version 1.0

workflow family_based_analysis {
  input {
    File input_vcf
    File input_vcf_index
    File pedigree
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    Array[String] chromosomes
    String callset_name
    String family_structure ## Trio or other
    
    File pLI_lookup
    File oe_lof_upper_lookup
    File oe_mis_upper_lookup
    File oe_syn_upper_lookup
    File clinvar_gene_desc
    File hgncSymbol_inheritance
    File gene_fullname_lookup
    File OMIM_entry_genesymbol

    # Runtime attributes
    String bcftools_docker = "quay.io/biocontainers/bcftools:1.13--h3a49de5_0"
    String SnpSift_docker = "quay.io/biocontainers/snpsift:4.3.1t--hdfd78af_3"
    String slivar_docker = "zihhuafang/slivar_modified:0.2.7"
    String vep_docker = "zihhuafang/ensembl_vep_loftee:v107"
    String report_docker = "zihhuafang/gp2_family_report:v1"
    String? gatk_docker_override
    String gatk_docker = select_first([gatk_docker_override, "broadinstitute/gatk:4.2.6.1"])
    String? gatk_path_override
    String gatk_path = select_first([gatk_path_override, "/gatk/gatk"])
    String? runtime_zones_override
    String runtime_zones = select_first([runtime_zones_override, "us-central1-a us-central1-b"])

    String? max_retries_override
    Int max_retries = select_first([max_retries_override, "1"])
    String? preemptible_tries_override
    Int preemptible_tries = select_first([preemptible_tries_override, "2"])
    Int? increase_disk_size
    Int additional_disk_space_gb = select_first([increase_disk_size, 20])
  }
    # Here we mesaure the size of the input vcf to decide the disk size of the VMs 
    Int disk_size = ceil((size(input_vcf, "GB")) * 1.5 + size(ref_fasta,"GB")) + additional_disk_space_gb


  parameter_meta {
    family_structure: "Sets the family structure type; accepted values are: `Trio` and `others`."
    input_vcf: "Set the path to your input vcf (one per family)"
    pedigree: "The Path of ped file (plink format) indicating relationship, sex and phenotype of family members."
    ref_fasta: "Sets the path to the reference."
    ref_fasta_index: "[Optional] Sets the path to the index of reference. If not provided, it will be inferred from the reference filename."
    chromosomes: "A list of the chromosomes for parallel processing for vep annotation."
    callset_name: "A string for the final output name, recommand to use the family ID."
  }

  scatter (chromosome in chromosomes) {
    call bcftools_norm_qc {
      input:
        input_vcf = input_vcf,
        input_vcf_index = input_vcf + ".tbi",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        chr = chromosome,
        docker = bcftools_docker,
        max_retries = max_retries,
        preemptible_tries = preemptible_tries,
        runtime_zones =runtime_zones,
        additional_disk_space_gb = additional_disk_space_gb
    }

    call vep {
      input:
        input_vcf = bcftools_norm_qc.output_vcf,
        input_vcf_index = bcftools_norm_qc.output_vcf_index,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        chr = chromosome,
        output_vcf = chromosome +"_vep.vcf.gz",
        docker = vep_docker,
        max_retries =max_retries,
        preemptible_tries = preemptible_tries,
        runtime_zones = runtime_zones,
        additional_disk_space_gb = additional_disk_space_gb
    }
  }

  call Gathervcf {
    input:
      input_vcfs_fofn = write_lines(vep.output_vcf),
      output_vcf_name = callset_name + "_vep.vcf.gz",
      disk_size = disk_size,
      docker = gatk_docker,
      gatk_path = gatk_path,
      max_retries = max_retries,
      preemptible_tries = preemptible_tries,
      runtime_zones = runtime_zones
  } 

  call slivar_seg {
    input:
      family_structure = family_structure,
      input_vcf = Gathervcf.output_vcf,
      input_vcf_index = Gathervcf.output_vcf_index,
      pedigree = pedigree,
      output_filename = callset_name + ".vcf",
      docker = slivar_docker,
      max_retries =max_retries,
      preemptible_tries = preemptible_tries,
      runtime_zones = runtime_zones
  }

  call Snpsift_count {
    input:
      input_vcf = slivar_seg.output_vcf,
      pedigree = pedigree,
      output_vcf_name = callset_name + ".vcf",
      docker = SnpSift_docker,
      max_retries =max_retries,
      preemptible_tries = preemptible_tries,
      runtime_zones = runtime_zones
  }

  call slivar_comphet {
    input:
      input_vcf = Snpsift_count.output_vcf,
      pedigree = pedigree,
      output_vcf_name = callset_name + ".vcf",
      docker = slivar_docker,
      max_retries =max_retries,
      preemptible_tries = preemptible_tries,
      runtime_zones = runtime_zones
  }

  call slivar_tsv {
    input:
      seg_vcf = slivar_seg.output_vcf,
      comphet_vcf = slivar_comphet.output_vcf,
      pedigree = pedigree,
      output_tsv_name = callset_name + ".tsv",
      pLI_lookup = pLI_lookup,
      oe_lof_upper_lookup = oe_lof_upper_lookup,
      oe_mis_upper_lookup = oe_mis_upper_lookup,
      oe_syn_upper_lookup = oe_syn_upper_lookup,
      clinvar_gene_desc =clinvar_gene_desc,
      hgncSymbol_inheritance = hgncSymbol_inheritance,
      gene_fullname_lookup = gene_fullname_lookup,
      OMIM_entry_genesymbol = OMIM_entry_genesymbol,
      docker = slivar_docker,
      max_retries =max_retries,
      preemptible_tries = preemptible_tries,
      runtime_zones = runtime_zones
  }

  call seg_var_report {
    input:
      input_tsv = slivar_tsv.output_tsv,
      docker = report_docker,
      max_retries =max_retries,
      preemptible_tries = preemptible_tries,
      runtime_zones = runtime_zones
  }

  call slivar_panel_filter {
    input:
      input_vcf = Gathervcf.output_vcf,
      input_vcf_index = Gathervcf.output_vcf_index,
      pedigree = pedigree,
      output_vcf_name = callset_name + ".vcf",
      docker = slivar_docker,
      max_retries =max_retries,
      preemptible_tries = preemptible_tries,
      runtime_zones = runtime_zones
  }

  call panel_tsv {
    input:
      input_vcf = slivar_panel_filter.output_vcf,
      pedigree = pedigree,
      output_tsv_name = callset_name + ".tsv",
      pLI_lookup = pLI_lookup,
      oe_lof_upper_lookup = oe_lof_upper_lookup,
      oe_mis_upper_lookup = oe_mis_upper_lookup,
      oe_syn_upper_lookup = oe_syn_upper_lookup,
      clinvar_gene_desc =clinvar_gene_desc,
      hgncSymbol_inheritance = hgncSymbol_inheritance,
      gene_fullname_lookup = gene_fullname_lookup,
      OMIM_entry_genesymbol = OMIM_entry_genesymbol,
      docker = slivar_docker,
      max_retries =max_retries,
      preemptible_tries = preemptible_tries,
      runtime_zones = runtime_zones
  }

  call panel_report {
    input:
      input_tsv = panel_tsv.output_tsv,
      output_prefix = callset_name,
      docker = report_docker,
      max_retries =max_retries,
      preemptible_tries = preemptible_tries,
      runtime_zones = runtime_zones
  }
  
  output {
    File segregation_vars_report = seg_var_report.output_excel
    File panel_var_report = panel_report.output_excel
    }
}

  
#
task bcftools_norm_qc {
  input {
    File input_vcf
    File input_vcf_index
    File ref_fasta
    File ref_fasta_index
    String chr
    Int threads = 5

    # Runtime parameters
    Int max_retries
    Int preemptible_tries
    String docker
    String runtime_zones
    Int additional_disk_space_gb
  }
    String output_filename = chr + "_norm.vcf.gz"
    Int disk_size = ceil((size(input_vcf, "GB"))* 2 + size(ref_fasta,"GB")) + additional_disk_space_gb
    
  command <<<
      # We subset variants that pass VQSR filter and by chromosome and perfom some basic variant QC    
   
       bcftools view -f PASS --regions ~{chr} -Ou ~{input_vcf} \
       | bcftools norm \
        -f ~{ref_fasta} \
        --threads ~{threads} \
        -m -any \
        -O u \
       | bcftools +fill-tags -O u -- -t AN,AC,AC_Hom,FORMAT/VAF \
       | bcftools view --trim-alt-alleles --min-ac 1 -O u \
       | bcftools filter -S . -i 'GQ >= 20 | DP >= 5 |(FMT/VAF >= 0.2 && FMT/VAF <= 0.8) | FMT/VAF <= 0.02 | FMT/VAF >= 0.98' -O z -o {output_filename}

        tabix -p vcf ~{output_filename}
    >>>

  runtime {
    docker: docker
    memory: "20 GB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: threads
    maxRetries: max_retries
    preemptible: preemptible_tries
    zones: runtime_zones 
  }

  output {
      File output_vcf = "~{output_filename}"
      File output_vcf_index = "~{output_filename}.tbi"
  }
}

task vep {
  input {
    File input_vcf
    File input_vcf_index
    File vep_cache
    File ref_fasta
    File ref_fasta_index
    String chr
    String output_vcf
    # Runtime parameters
    Int additional_disk_space_gb
    Int max_retries
    Int preemptible_tries
    String docker
    String runtime_zones
  }
    
    Int disk_size = ceil((size(input_vcf, "GB"))* 3 + size(ref_fasta,"GB") + size(vep_cache, "GB")*2) + additional_disk_space_gb
    String output_filename = chr + "_norm.vcf.gz"

  command <<<

  tar xzf ~{vep_cache}
  VEP_CACHE=$(basename ~{vep_cache} .tar)

    vep \
      -i ~{input_vcf} \
      --ASSEMBLY GRCh38 \
      --buffer_size 50000 \
      --fasta ~{ref_fasta} \
      --sift b \
      --polyphen b \
      --ccds \
      --hgvs \
      --symbol \
      --numbers \
      --domains \
      --regulatory \
      --canonical \
      --protein \
      --biotype \
      --af \
      --af_1kg \
      --af_gnomade \
      --max_af \
      --pubmed \
      --uniprot \
      --mane \
      --tsl \
      --appris \
      --variant_class \
      --gene_phenotype \
      --mirna \
      --var_synonyms \
      --check_existing \
      --nearest symbol \
      --terms SO \
      --check_existing \
      --clin_sig_allele 1 \
      --force_overwrite \
      --cache \
      --pick \
      --pick_order biotype,canonical,appris,tsl,ccds,rank,length,mane \
      --offline \
      --dir_cache ${VEP_CACHE} \
      --dir_plugins ${VEP_CACHE}/Plugins/ \
      --plugin UTRannotator,${VEP_CACHE}/vep_annotation/uORF_5UTR_GRCh38_PUBLIC.txt \
      --plugin SpliceAI,snv=${VEP_CACHE}/vep_annotation/spliceai_scores.masked.indel.hg38.vcf.gz,indel=vep_cache/vep_annotation/spliceai_scores.masked.indel.hg38.vcf.gz \
      --plugin SpliceRegion,Extended \
      --plugin TSSDistance \
      --plugin NMD \
      --plugin Downstream \
      --plugin DisGeNET,file=${VEP_CACHE}/vep_annotation/all_variant_disease_pmid_associations_final.tsv.gz,disease=1 \
      --plugin CADD,
      --plugin dbNSFP,${VEP_CACHE}/vep_annotation/dbNSFP4.3a_grch38.gz,Ensembl_transcriptid,MetaRNN_score,MPC_score,LRT_score,GERP++_RS,FATHMM_score,fathmm-MKL_coding_score,DANN_score,REVEL_score,PrimateAI_score,MutPred_score,GTEx_V8_gene,GTEx_V8_tissue,Geuvadis_eQTL_target_gene,gnomAD_genomes_controls_and_biobanks_nhomalt \
      --plugin LoF,loftee_path:${VEP_CACHE}/Plugins/loftee_hg38,\
        human_ancestor_fa:${VEP_CACHE}/Plugins/loftee_hg38/human_ancestor.fa.gz,\
conservation_file:${VEP_CACHE}/Plugins/loftee_hg38/loftee.sql,\
gerp_bigwig:${VEP_CACHE}/Plugins/loftee_hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
      --custom ${VEP_CACHE}/vep_annotation/clinvar_20220730.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNDN \
      -o ~{output_vcf} \
      --compress_output bgzip \
      --vcf \
      --stats_text
  >>>

  runtime {
    docker: docker
    memory: "10 GB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: 1
    maxRetries: max_retries
    preemptible: preemptible_tries
    zones: runtime_zones 
  }

  output {
      File output_vcf = "~{output_filename}"
      File output_vcf_index = "~{output_filename}.tbi"
  }
}

task Gathervcf {
  input {
    File input_vcfs_fofn
    String output_vcf_name
    String gatk_path

    String docker
    String runtime_zones
    Int disk_size
    Int max_retries
    Int preemptible_tries
  }
    

  command <<<
    set -e

    # Now using NIO to localize the vcfs but the input file must have a ".list" extension
    mv ~{input_vcfs_fofn} inputs.list

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    ${gatk_path} --java-options "-Xmx15g -Xms15g" \
    GatherVcfsCloud \
    --ignore-safety-checks \
    --gather-type BLOCK \
    --input inputs.list \
    --output ~{output_vcf_name}

    ${gatk_path} --java-options "-Xmx6g -Xms6g" \
    IndexFeatureFile \
    --input ~{output_vcf_name}
  >>>
  runtime {
    docker: docker
    memory: "20 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    maxRetries: max_retries
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}

task slivar_seg {
  input {
    File input_vcf
    File input_vcf_index 
    File pedigree
    String output_filename
    String family_structure

    String docker
    Int max_retries
    Int preemptible_tries
    String runtime_zones
  }

    Int disk_size = ceil(size(input_vcf, "GB") * 3)

  command <<<
    set -e
        
    if [[ "~{family_structure}" == "Trio" ]]; then
        slivar expr --vcf ~{input_vcf} \
                    --ped ~{pedigree} \
                    -o ~{output_filename} \
                    --js /opt/slivar/slivar-functions.js \
                    --pass-only \
                    --info "variant.ALT[0] != '*' && variant.call_rate >= 0.95" \
                    --trio "denovo:(variant.CHROM != 'X' || variant.CHROM != 'chrX') && denovo(kid, mom, dad) && INFO.gnomADg_Popmax_AF < 0.01" \
                    --trio "x_denovo:x_denovo(kid, mom, dad) && INFO.gnomADg_Popmax_AF <= 0.01" \
                    --trio "HOM:(variant.CHROM != 'X' || variant.CHROM != 'chrX') && recessive(kid, mom, dad) && INFO.gnomADg_nhomalt < 10" \
                    --trio "X_HOM:x_recessive(kid, mom, dad) && INFO.gnomADg_nhomalt <= 10" \
                    --trio "comphet_side:comphet_side(kid, mom, dad) && INFO.gnomADg_nhomalt <= 10" \
                    --family-expr "dominant:(variant.CHROM != 'X' || variant.CHROM != 'chrX') && fam.every(segregating_dominant) && !fam.every(segregating_denovo) && INFO.gnomADg_Popmax_AF <= 0.01" \
                    --family-expr "x_dominant:(variant.CHROM == 'X' || variant.CHROM == 'chrX') && fam.every(segregating_dominant_x) && !fam.every(segregating_denovo_x) && INFO.gnomADg_Popmax_AF <= 0.01" 
    else
        slivar expr --vcf ~{input_vcf} \
                    --ped ~{pedigree} \
                    -o ~{output_filename} \
                    --js /opt/slivar/slivar-functions.js \
                    --pass-only \
                    --info "variant.ALT[0] != '*' && variant.call_rate >= 0.95" \
                    --family-expr "HOM:(variant.CHROM != 'X' || variant.CHROM != 'chrX') && fam.every(segregating_recessive) && INFO.gnomADg_nhomalt <= 10" \
                    --family-expr "X_HOM:(variant.CHROM == 'X' || variant.CHROM == 'chrX') && fam.every(segregating_recessive_x) && INFO.gnomADg_nhomalt <= 10" \
                    --family-expr "denovo:(variant.CHROM != 'X' || variant.CHROM != 'chrX') && fam.every(segregating_denovo) && INFO.gnomADg_Popmax_AF <= 0.01" \
                    --family-expr "x_denovo:(variant.CHROM == 'X' || variant.CHROM == 'chrX') && fam.every(segregating_denovo_x) && INFO.gnomADg_Popmax_AF <= 0.01" \
                    --family-expr "dominant:(variant.CHROM != 'X' || variant.CHROM != 'chrX') && fam.every(segregating_dominant) && !fam.every(segregating_denovo) && INFO.gnomADg_Popmax_AF <= 0.01" \
                    --family-expr "x_dominant:(variant.CHROM == 'X' || variant.CHROM == 'chrX') && fam.every(segregating_dominant_x) && !fam.every(segregating_denovo_x) && INFO.gnomADg_Popmax_AF <= 0.01" \
                    --family-expr "comphet_side:fam.every(function(s) {{return (s.het || s.hom_ref) && hq1(s) }}) && fam.some(function(s) {return s.het && s.affected}) && INFO.gnomADg_nhomalt <= 10"
    fi

  >>>
  runtime {
    docker: docker
    memory: "20 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    maxRetries: max_retries
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {
    File output_vcf = "~{output_filename}"
  }
}

task Snpsift_count {
  input {
    File input_vcf
    File pedigree
    String output_vcf_name

    String docker
    Int max_retries
    Int preemptible_tries
    String runtime_zones
  }
    Int disk_size = ceil(size(input_vcf, "GB") * 2.5)

  command <<<
    set -e
    mv ~{pedigree} fam.tfam
    SnpSift caseControl -tfam fam.tfam ~{input_vcf} > {output_vcf_name}

  >>>
  runtime {
    docker: docker
    memory: "20 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    maxRetries: max_retries
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {
    File output_vcf = "~{output_vcf_name}"
  }
}

task slivar_comphet {
  input {
    File input_vcf
    File pedigree
    String output_vcf_name

    String docker
    Int max_retries
    Int preemptible_tries
    String runtime_zones
  }
    Int disk_size = ceil(size(input_vcf, "GB") * 2.5)

  command <<<
    set -e

    slivar compound-hets \
           -v {input_vcf} 
           --allow-non-trios --sample-field comphet_side --sample-field denovo -p ~{pedigree} \
           -o {output_vcf_name}
           
  >>>
  runtime {
    docker: docker
    memory: "3.75 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    maxRetries: max_retries
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {
    File output_vcf = "~{output_vcf_name}"
  }
}

# Creating a spreadsheet from a filtered VCF with variants fit differnet mode of inheritance
task slivar_tsv {
  input {
    File seg_vcf
    File comphet_vcf
    File pedigree
    String output_tsv_name

    String docker
    Int max_retries
    Int preemptible_tries
    String runtime_zones

    File pLI_lookup
    File oe_lof_upper_lookup
    File oe_mis_upper_lookup
    File oe_syn_upper_lookup
    File clinvar_gene_desc
    File hgncSymbol_inheritance
    File gene_fullname_lookup
    File OMIM_entry_genesymbol
  }
    Int disk_size = ceil(size(seg_vcf, "GB") * 2.5)
    
    Array[String] info_field = ["gnomADg_AF", "gnomADg_Popmax_AF", "gnomADg_nhomalt", "gnomADg_AC", "gnomADg_Popmax_AF_filter", "TOPMed8_AF", "CADD_PHRED", "Cases", "Controls"]

    Array[String] csq_columns= [
    "Existing_variation",
    "gnomADe_AF",
    "MAX_AF_POPS",
    "LoF",
    "LoF_filter",
    "LoF_flags",
    "LoF_info",
    "MetaRNN_score",
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",
    "SpliceAI_pred_SYMBOL",
    "SpliceRegion",
    "existing_InFrame_oORFs",
    "existing_OutOfFrame_oORFs",
    "existing_uORFs",
    "five_prime_UTR_variant_annotation",
    "five_prime_UTR_variant_consequence",
    "Gene",
    "IMPACT",
    "BIOTYPE",
    "STRAND",
    "CANONICAL",
    "MANE_SELECT",
    "TSL",
    "NEAREST",
    "EXON",
    "Codons",
    "Amino_acids",
    "HGVSc",
    "HGVSp",
    "ClinVar_CLNSIG",
    "ClinVar_CLNDN",
    "VAR_SYNONYMS",
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "TSSDistance",
    "miRNA",
    "TRANSCRIPTION_FACTORS",
    "NMD",
    "DisGeNET",
    "DownstreamProtein",
    "ProteinLengthChange"
    ]

  command <<<
    set -e

    # Variants that fit (x-) domoinant, (x-denovo) and (x-) HOM (homozygotes) 
    slivar tsv \
        --info-field ~{sep=" --info-field " info_field} \
        -s denovo \
        -s x_denovo \
        -s HOM \
        -s X_HOM \
        -s dominant \
        -s x_dominant \
        -c CSQ \
        --csq-column ~{sep=" --csq-column " csq_columns} \
        -g ~{pLI_lookup}
        -g ~{oe_lof_upper_lookup} \
        -g ~{oe_mis_upper_lookup} \
        -g ~{oe_syn_upper_lookup} \
        -g ~{clinvar_gene_desc} \
        -g ~{hgncSymbol_inheritance} \
        -g ~{gene_fullname_lookup} \
        -g ~{OMIM_entry_genesymbol} \
        -p ~{pedigree} \
        ~{seg_vcf} > seg.tsv

    # Compound heterzygotes
    slivar tsv \
        -s slivar_comphet \
        --info-field ~{sep=" --info-field " info_field} \
        -c CSQ \
        --csq-column ~{sep=" --csq-column " csq_columns} \
        -g ~{pLI_lookup}
        -g ~{oe_lof_upper_lookup} \
        -g ~{oe_mis_upper_lookup} \
        -g ~{oe_syn_upper_lookup} \
        -g ~{clinvar_gene_desc} \
        -g ~{hgncSymbol_inheritance} \
        -g ~{gene_fullname_lookup} \
        -g ~{OMIM_entry_genesymbol} \
        -p ~{pedigree} \
        ~{comphet_vcf} | { grep -v ^# || true; } >> comphet.tsv

    # We merge two tsv together to make the same spreadsheet
    cat seg.tsv comphet.tsv | sed -r -e 's/slivar_comphet_[0-9]+/comphet/g' \
         | sed '1 s/gene_description_1/gnomAD_pLI/;s/gene_description_2/gnomAD_oe_lof_upper/;s/gene_description_3/gnomAD_oe_mis_upper/;s/gene_description_4/gnomAD_oe_syn_upper/;s/gene_description_5/ClinVar_gene_description/;s/gene_description_6/MOI/;s/gene_description_7/Gene_Fullname/;s/gene_description_8/OMIM/;' > ~{output_tsv_name}
  

  >>>
  runtime {
    docker: docker
    memory: "3.75 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    maxRetries: max_retries
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {
    File output_tsv = "~{output_tsv_name}"
  }
}

task seg_var_report {
  input {
    File input_tsv
    String output_prefix
    
    String docker
    Int max_retries
    Int preemptible_tries
    String runtime_zones
  }
    Int disk_size = ceil(size(input_tsv, "GB") * 1.5)


  command <<<
    set -e
    
    python /opt/bin/report_format.py ~{input_tsv} ~{output_prefix}

  >>>
  runtime {
    docker: docker
    memory: "3.75 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    maxRetries: max_retries
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {
    File output_excel = "~{output_prefix}.xlsx"
  }
}

task slivar_panel_filter {
  input {
    File input_vcf
    File input_vcf_index 
    File pedigree
    File input_bed
    String output_vcf_name

    String docker
    Int max_retries
    Int preemptible_tries
    String runtime_zones
  }
    Int disk_size = ceil(size(input_vcf, "GB") * 2.5)

  command <<<
    set -e

    slivar expr --vcf ~{input_vcf} \
                    -o ~{output_vcf_name} \
                    --ped ~{pedigree} \
                    --js /opt/slivar/slivar-functions.js \
                    --pass-only \
                    --info 'INFO.gnomad_AF < 0.05 && variant.ALT[0] != "*" && variant.call_rate >= 0.95' \
                    --region ~{input_bed} \
                    --family-expr "HOM:fam.every(function(s) {return s.hom_alt && s.affected && hq1(s)})" \
                    --family-expr "HET:fam.every(function(s) {return s.het && s.affected && hq1(s)})"

  >>>
  runtime {
    docker: docker
    memory: "10 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    maxRetries: max_retries
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {
    File output_vcf = "~{output_vcf_name}"
  }
}

task panel_tsv {
  input {
    File input_vcf
    File pedigree
    String output_tsv_name

    String docker
    Int max_retries
    Int preemptible_tries
    String runtime_zones
    
    File pLI_lookup
    File oe_lof_upper_lookup
    File oe_mis_upper_lookup
    File oe_syn_upper_lookup
    File clinvar_gene_desc
    File hgncSymbol_inheritance
    File gene_fullname_lookup
    File OMIM_entry_genesymbol
  }  
    Int disk_size = ceil(size(input_vcf, "GB") * 2.5)

    Array[String] info_field =[
    "gnomADg_AF",
    "gnomADg_Popmax_AF",
    "gnomADg_nhomalt",
    "gnomADg_AC",
    "gnomADg_Popmax_AF_filter",
    "TOPMed8_AF",
    "CADD_PHRED"
    ]

    Array[String] csq_columns= [
    "Existing_variation",
    "gnomADe_AF",
    "MAX_AF_POPS",
    "LoF",
    "LoF_filter",
    "LoF_flags",
    "LoF_info",
    "MetaRNN_score",
    "SpliceAI_pred_DS_AG",
    "SpliceAI_pred_DS_AL",
    "SpliceAI_pred_DS_DG",
    "SpliceAI_pred_DS_DL",
    "SpliceAI_pred_SYMBOL",
    "SpliceRegion",
    "existing_InFrame_oORFs",
    "existing_OutOfFrame_oORFs",
    "existing_uORFs",
    "five_prime_UTR_variant_annotation",
    "five_prime_UTR_variant_consequence",
    "Gene",
    "IMPACT",
    "BIOTYPE",
    "STRAND",
    "CANONICAL",
    "MANE_SELECT",
    "TSL",
    "NEAREST",
    "EXON",
    "Codons",
    "Amino_acids",
    "HGVSc",
    "HGVSp",
    "ClinVar_CLNSIG",
    "ClinVar_CLNDN",
    "VAR_SYNONYMS",
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "TSSDistance",
    "miRNA",
    "TRANSCRIPTION_FACTORS",
    "NMD",
    "DisGeNET",
    "DownstreamProtein",
    "ProteinLengthChange"
    ]
    
  command <<<
    set -e

    # Variants that are homozygous or heterozygous in gene panels
    slivar tsv \
        --info-field ~{sep=" --info-field " info_field} \
        -s HOM \
        -s HET \
        -c CSQ \
        --csq-column ~{sep=" --csq-column " csq_columns} \
        -g ~{pLI_lookup}
        -g ~{oe_lof_upper_lookup} \
        -g ~{oe_mis_upper_lookup} \
        -g ~{oe_syn_upper_lookup} \
        -g ~{clinvar_gene_desc} \
        -g ~{hgncSymbol_inheritance} \
        -g ~{gene_fullname_lookup} \
        -g ~{OMIM_entry_genesymbol} \
        -p ~{pedigree} \
        ~{input_vcf} > tmp.tsv


    # Modify the header for slivar gene description 
    cat tmp.tsv \
         | sed '1 s/gene_description_1/gnomAD_pLI/;s/gene_description_2/gnomAD_oe_lof_upper/;s/gene_description_3/gnomAD_oe_mis_upper/;s/gene_description_4/gnomAD_oe_syn_upper/;s/gene_description_5/ClinVar_gene_description/;s/gene_description_6/MOI/;s/gene_description_7/Gene_Fullname/;s/gene_description_8/OMIM/;' > ~{output_tsv_name}
  

  >>>
  runtime {
    docker: docker
    memory: "3.75 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    maxRetries: max_retries
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {
    File output_tsv = "~{output_tsv_name}"
  }
}

task panel_report {
  input {
    File input_tsv
    File core_panel_gene_list
    File extend_panel_gene_list
    String output_prefix
    
    String docker
    Int max_retries
    Int preemptible_tries
    String runtime_zones
  }
    Int disk_size = ceil(size(input_tsv, "GB") * 1.5)


  command <<<
    set -e
    
    python /opt/bin/report_panel.py ~{input_tsv} ~{core_panel_gene_list} ~{extend_panel_gene_list} ~{output_prefix}

  >>>
  runtime {
    docker: docker
    memory: "3.75 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    maxRetries: max_retries
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {
    File output_excel = "~{output_prefix}.xlsx"
  }
}
