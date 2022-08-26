version 1.0

# The WDL performs joint-genotyping of sv discovered from parliament2 with muCNV
# in human whole-genome sequencing (WGS)
#
# Input: individual *.survivor.qual.vcf from the output of parliament2
# Output: a joint-called SV set

# WORKFLOW DEFINITION
workflow muCNV {
  input {
    Array[File] input_vcfs
    Array[File] input_crams
    Int num_crams = length(input_crams)
    String ref_gc = "/muCNV/resource/GRCh38.v3.gc"
    File refFasta
    File refFasta_index
    File refFasta_dict
    File ref_cache

    # Runtime paramters
    String docker_survivor = "szarate/survivor:1d1d33b"
    String docker_muCNV = "avalillarionova/mucnv:1.0.0"
    String docker_bcftools = "quay.io/biocontainers/bcftools:1.13--h3a49de5_0"
    String runtime_zones 
    Int max_retries = 1
    Int preemptible_tries = 2
    
    Int? increase_disk_size
    Int additional_disk_space_gb = select_first([increase_disk_size, 20])
    Int? small_disk_override
    Int small_disk = select_first([small_disk_override, "100"])
    Int? medium_disk_override
    Int medium_disk = select_first([medium_disk_override, "200"])
    Int? large_disk_override
    Int large_disk = select_first([large_disk_override, "300"])
    Int? huge_disk_override
    Int huge_disk = select_first([huge_disk_override, "400"])    
  }
  
  call prep_chr {
    input:
      refFasta_index = refFasta_index,
      runtime_zones = runtime_zones
    }
  
  Array[String] chrs = read_lines(prep_chr.chrs)
  
  scatter (idx in range(length(input_vcfs))) { 
 
  # keep only PASS SV from parliament output
    call filtersv {
      input:
        input_vcf = input_vcfs[idx],
        sample=basename(input_vcfs[idx], ".survivor.qual.vcf"),
        docker_bcftools = docker_survivor,
        runtime_zones = runtime_zones,
        preemptible_tries = preemptible_tries,
        additional_disk_space_gb = additional_disk_space_gb
    }
  }
  
  scatter (chr in chrs) {
    call vcf_fofn {
      input:
        fofn = write_lines(flatten(filtersv.filtered_vcfs)),
        chr = chr        
    }  
  }
  
  scatter (fofn in vcf_fofn.file_list) {
    # Merge individual filtered sv
    call mergeSV as mergeSVchrom {
      input:
        input_vcfs = read_lines(fofn),
        chr = sub(basename(fofn, ".list"), "_fofn",""),
        docker_survivor = docker_survivor,
        runtime_zones =runtime_zones,
        max_retries = max_retries,
        preemptible_tries = preemptible_tries,
        additional_disk_space_gb = large_disk
    }
  }
  
  call vcfconcat {
     input:
       input_vcfs = mergeSVchrom.survivorSortedVCF,
       docker = docker_survivor,
       runtime_zones = runtime_zones,
       preemptible_tries = preemptible_tries,
       max_retries = max_retries,
       vm_disk_size = huge_disk
     }
  
  # Prepare vcf interval from merged vcf for muCNV
    
  call prep_interval {
    input:
      input_vcf = vcfconcat.output_vcf,
      docker_muCNV = docker_muCNV,
      runtime_zones = runtime_zones,
      max_retries = max_retries,
      preemptible_tries = preemptible_tries,
      additional_disk_space_gb = additional_disk_space_gb
    }
    
  # muCNV pileup per sample 
  scatter (idx in range(length(input_crams))) { 

    call pileup {
      input:
        input_cram=input_crams[idx],
        input_cram_index=input_crams[idx]+ ".crai",
        sample=basename(input_crams[idx], ".cram"),
        vcf_interval = prep_interval.out_interval,
        refFasta = refFasta,
        refFasta_index = refFasta_index, 
        refFasta_dict = refFasta_dict,
        ref_cache = ref_cache,
        ref_gc = ref_gc,
        docker_muCNV = docker_muCNV,
        runtime_zones = runtime_zones,
        max_retries = max_retries,
        preemptible_tries = preemptible_tries,
        additional_disk_space_gb = additional_disk_space_gb
    }
  }

  # For small callsets (fewer than 500 samples) we can do joint-genotpying directly.
  # For anything larger, we need to perform muCNV merge first.
  
  Boolean is_small_callset = num_crams <= 500
  
  if (is_small_callset) {
    call final_jg {
      input:
        vcf_interval = prep_interval.out_interval,
        input_pileup = pileup.pileup_out,
        input_idx = pileup.idx_out,
        input_var = pileup.var_out,
        ref_gc = ref_gc, 
        docker_muCNV = docker_muCNV,
        runtime_zones = runtime_zones,
        max_retries = max_retries,
        preemptible_tries = preemptible_tries,
        additional_disk_space_gb = additional_disk_space_gb
    }
  }
    # joint-genotyping sv
  if (!is_small_callset) {
    
    call muCNV_merge{
      input:
        vcf_interval = prep_interval.out_interval,
        input_pileup = pileup.pileup_out,
        input_idx = pileup.idx_out,
        input_var = pileup.var_out,
        ref_gc = ref_gc,
        docker_muCNV = docker_muCNV,
        runtime_zones = runtime_zones,
        max_retries = max_retries,
        preemptible_tries = preemptible_tries,
        additional_disk_space_gb = additional_disk_space_gb
    }

    scatter (idx in range(length(muCNV_merge.merge_pileup))) {
      
      call merge_jg {
        input:
          vcf_interval = prep_interval.out_interval,
          input_pileup = muCNV_merge.merge_pileup[idx],
          input_var = muCNV_merge.merge_var[idx],
          input_idx = muCNV_merge.merge_idx[idx],
          ref_gc = ref_gc,
          chromosome = sub(basename(muCNV_merge.merge_pileup[idx], ".pileup"),"GP2.",""),
          docker_muCNV = docker_muCNV,
          runtime_zones = runtime_zones,
          max_retries = max_retries,
          preemptible_tries = preemptible_tries,
          additional_disk_space_gb = additional_disk_space_gb
      }
    }
  }
  # Get the VCFs from either code path
  #Array[File?] output_vcf_files = if defined(final_jg.outvcf) then [final_jg.outvcf] else merge_jg.merge_outvcf
  #Array[File?] output_vcf_index_files = if defined(final_jg.outvcf_index) then [final_jg.outvcf_index] else merge_jg.merge_outvcf_index

  output {
    # output the output of Jasmine
    File nongeno_vcf = vcfconcat.output_vcf

    # outputs from the small callset or large callset through the wdl
    #Array[File?] output_vcfs = select_all(output_vcf_files)
    #Array[File?] output_vcf_indices = select_all(output_vcf_index_files)
    
    #output from small callset
    File? small_callset_vcf = final_jg.outvcf
    
    #output from large callset
    Array[File]? large_callset_vcf = merge_jg.merge_outvcf
  }
}

#__________________________________________________________________________________________________________

#TASK DEFINITIONS
# prepare chromosome name
task prep_chr {
  input {
    # Command parameters
    File refFasta_index
    # Runtime parameters
    String runtime_zones 
  }

    command <<<
        
        ## we sort the order cuz glob in muCNV_merge is lexicographic
        cat ~{refFasta_index} | head -n 22 | awk '{print $1}' | sort > chrs.txt
  
    >>>

    runtime {
        docker : "ubuntu:20.04"
        memory: "5 GB"
        cpu: 1
        zones: runtime_zones
    }

    output {
        File chrs = "chrs.txt"
    }
}


# Keep only PASS sv from parliament output
task filtersv {
  input {
    # Command parameters
    File input_vcf
    String sample
    
    # Runtime parameters
    String docker_bcftools
    String runtime_zones 
    Int preemptible_tries
    Int additional_disk_space_gb
  }
    String output_filename = sample + ".vcf"
    Int disk_size = ceil(size(input_vcf, "GB") * 2 ) + additional_disk_space_gb

    command <<<
        
    for i in {1..22}
    do
        vcftools --vcf ~{input_vcf} --remove-filtered-all --chr "chr${i}" --recode --recode-INFO-all --stdout > "chr${i}_~{output_filename}"
    done
    
    >>>

    runtime {
        docker : docker_bcftools
        memory: "5 GB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        zones: runtime_zones
        preemptible: preemptible_tries
    }

    output {
       Array[File] filtered_vcfs = glob("chr*_~{output_filename}")
    }
}

task vcf_fofn {
  input {
    File fofn
    String chr
  }
    command <<<
        
    echo "~{chr}"    
    
    grep '~{chr}_' ~{fofn} > ~{chr}_fofn.list
  
    >>>

    runtime {
        docker : "ubuntu:20.04"
        memory: "3.75 GB"
        cpu: 1
    }
    
    output {
        File file_list = "~{chr}_fofn.list"
    }
}


# Merge individual filtered sv from parliament output
task mergeSV {
  input {
    # Command parameters
    Array[File] input_vcfs
    String chr
    # Runtime parameters
    Int additional_disk_space_gb
    Int max_retries
    Int preemptible_tries
    String docker_survivor
    String runtime_zones 
  }
    Int disk_size = ceil(size(input_vcfs, "GB") * 4) + additional_disk_space_gb
    String output_filename = chr + ".vcf"

    command <<<
        echo "~{disk_size}"
        echo "~{chr}" 
       
        echo "~{sep='\n' input_vcfs}" > survivor_inputs
        
        #We consider maximum allowed distance of 200 bp from muCNV paper
        #Report sv call present in at least one sample and agree on type and strand
        SURVIVOR merge survivor_inputs 200 1 1 0 0 10 tmp.vcf

        vcf-sort -c < tmp.vcf > "~{output_filename}"

    >>>

    runtime {
        docker : docker_survivor
        memory: "30 GB"
        disks: "local-disk " + disk_size + " HDD"
        zones: runtime_zones
        maxRetries: max_retries
        preemptible: preemptible_tries
    }

    output {
        File survivorSortedVCF = "~{output_filename}"
    }
}

task vcfconcat {
  input {
    # Command parameters
    Array[File] input_vcfs
    String output_filename = "gp2.vcf"

    # Runtime parameters
    String docker
    String runtime_zones
    Int vm_disk_size
    Int max_retries
    Int preemptible_tries
  }
  
    command <<<
    
    echo "~{vm_disk_size}"
    
    echo "~{sep='\n' input_vcfs}" |  awk -F '/' '{print $NF,$0}' | sort -V | cut -f2- -d" " > input_vcf.list
    
    vcf-concat -f input_vcf.list > ~{output_filename}

    >>>

  runtime {
    docker: docker
    disks: "local-disk " + vm_disk_size + " HDD"
    zones: runtime_zones
    preemptible: preemptible_tries
    memory: "20 GB"
    cpu: "1"
  }
  output {
    File input_list = "input_vcf.list"
    File output_vcf = "~{output_filename}"
  }
}
  
# Prepare interval input for mnCNV
task prep_interval {
  input {
    # Command parameters
    File input_vcf

    # Runtime parameters
    Int additional_disk_space_gb
    Int max_retries
    Int preemptible_tries
    String docker_muCNV
    String runtime_zones 
    String output_filename = "candidate_sv.interavls"
  }
    Int disk_size = ceil(size(input_vcf, "GB") * 2 ) + additional_disk_space_gb

    command <<<
        muCNV vcf2int -v ~{input_vcf} -i ~{output_filename}
    >>>

    runtime {
        docker : docker_muCNV
        memory: "8 GB"
        cpu: 1
        disks: "local-disk " + disk_size + " HDD"
        zones: runtime_zones
        maxRetries: max_retries
        preemptible: preemptible_tries
    }

    output {
        File out_interval = "~{output_filename}"
    }
}

# Generate pileup from cram for individual samples
task pileup {
  input {
    # Command parameters
    File input_cram
    File input_cram_index
    File vcf_interval
    String sample
    String ref_gc 
    File refFasta
    File refFasta_index
    File refFasta_dict
    File ref_cache
    # Runtime parameters
    String docker_muCNV
    Int max_retries
    Int preemptible_tries
    Int additional_disk_space_gb
    String runtime_zones 
  }

    Int disk_size = ceil(size(input_cram, "GB") * 2 + size(refFasta, "GB")) + additional_disk_space_gb

  command <<<
    
    echo "Untarring reference cache source"
    tar xzf ~{ref_cache}

    echo "Specifying paths to reference caches"
    export REF_PATH="$(pwd)/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
    export REF_CACHE="$(pwd)/ref/cache/%2s/%2s/%s"

    mkdir "~{sample}"
    
    muCNV pileup -s ~{sample} -V ~{vcf_interval} -b ~{input_cram} -f ~{ref_gc} -o "~{sample}"
    
    mv "~{sample}/~{sample}.pileup" "~{sample}.pileup"
    mv "~{sample}/~{sample}.idx" "~{sample}.idx"
    mv "~{sample}/~{sample}.var" "~{sample}.var"

    >>>

  output {
    File pileup_out = "~{sample}.pileup"
    File idx_out = "~{sample}.idx"
    File var_out = "~{sample}.var"
    #Array[File] pileup_out = glob("~{sample}.*")
  }

  runtime {
    docker: docker_muCNV
    memory: "10 GB"
    disks: "local-disk " + disk_size + " HDD"
    zones: runtime_zones
    maxRetries: max_retries
    preemptible: preemptible_tries
  }
}

## Merge pileup per chr when n of samples > 500
task muCNV_merge {
  input {
    # Command parameters
    Array[File] input_pileup
    Array[File] input_idx
    Array[File] input_var
    File vcf_interval
    String ref_gc 
    String output_prefix = "GP2"
    # Runtime parameters
    String docker_muCNV
    Int max_retries
    Int preemptible_tries
    Int additional_disk_space_gb
    String runtime_zones 
  }
    
    Int disk_size = ceil(size(input_pileup, "GB")*3) + additional_disk_space_gb

  command <<<
        
    echo "~{sep='\n' input_pileup}" | awk -F'.' '{print $1}' > pileup.list

    muCNV merge -V ~{vcf_interval} -i pileup.list -f ~{ref_gc} -o ~{output_prefix}

    >>>

  output {
    File pileup_list = "pileup.list"
    Array[File] merge_pileup = glob("~{output_prefix}.*.pileup")
    Array[File] merge_var = glob("~{output_prefix}.*.var")
    Array[File] merge_idx = glob("~{output_prefix}.*.idx")
  }

  runtime {
    docker: docker_muCNV
    memory: "20 GB"
    disks: "local-disk " + disk_size + " HDD"
    zones: runtime_zones
    maxRetries: max_retries
    preemptible: preemptible_tries
  }
}

task final_jg {
  input {
    # Command parameters
    File vcf_interval
    Array[File] input_pileup
    Array[File] input_idx
    Array[File] input_var
    String ref_gc 
    String output_filename = "gp2_sv.vcf"

    # Runtime parameters
    String docker_muCNV
    Int max_retries
    Int preemptible_tries
    Int additional_disk_space_gb
    String runtime_zones 

  }

    Int disk_size = ceil(size(input_pileup, "GB")) + additional_disk_space_gb

  command <<<
  
    mkdir -p pileup
    mv ~{input_pileup} ~{input_idx} ~{input_var} pileup/
    ls -a pileup/*.pileup | cat | awk -F'.' '{print $1}' > pileup.list

    muCNV genotype -V ~{vcf_interval} -i pileup.list -f ~{ref_gc} -o ~{output_filename}
    >>>

  output {
    File pileup_list = "pileup.list"
    File outvcf = "~{output_filename}"
  }

  runtime {
    docker: docker_muCNV
    memory: "20 GB"
    cpu: 1
    disks: "local-disk " + disk_size + " HDD"
    zones: runtime_zones
    maxRetries: max_retries
    preemptible: preemptible_tries
  }
}


task merge_jg {
  input {
    # Command parameters
    File vcf_interval
    File input_pileup
    File input_idx
    File input_var
    String ref_gc 
    String output_filename = chromosome + "_gp2_sv.vcf"
    String chromosome
    # Runtime parameters
    String docker_muCNV
    Int max_retries
    Int preemptible_tries
    Int additional_disk_space_gb
    String runtime_zones 

  }
    Int disk_size = ceil(size(input_pileup, "GB")) + additional_disk_space_gb
    String chr_string = sub(chromosome, "^chr", "")

  command <<<
    # Check the chromosome and the pileup files match
    echo "~{chromosome}"
    echo "~{input_pileup}"
    
    mkdir -p pileup
    mv ~{input_pileup} ~{input_idx} ~{input_var} pileup/
    
    echo "pileup/GP2" > pileup.list

    muCNV genotype -V ~{vcf_interval} -i pileup.list -c ~{chr_string} -f ~{ref_gc} -o ~{output_filename}

    >>>

  output {
    File merge_outvcf = "~{output_filename}"
    File merge_pileup_list = "pileup.list"
  }

  runtime {
    docker: docker_muCNV
    memory: "20 GB"
    disks: "local-disk " + disk_size + " HDD"
    zones: runtime_zones
    maxRetries: max_retries
    preemptible: preemptible_tries
  }
}