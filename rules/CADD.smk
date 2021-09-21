#!/usr/bin/env python
import glob
from pathlib import Path

parentdir = Path(srcdir("")).parents[0]
configfile: "configs/config_GRCh38_v1.6_noanno.yml"

rule prepare:
    input: 'CADD_input/{chromosome}.vcf'
    output: temp('CADD_tmp/{chromosome}.prepared.vcf')
    conda: '../envs/CADD.yml'
    shell:
        '''
        cat {input} \
        | python CADD_src_scripts/VCF2vepVCF.py \
        | sort -k1,1 -k2,2n -k4,4 -k5,5 \
        | uniq > {output}
        '''        

rule prescore:
    input: rules.prepare.output
    output:
        novel=temp('CADD_tmp/{chromosome}.novel.vcf'),
        prescored=temp('CADD_tmp/{chromosome}.pre.tsv')
    conda: '../envs/CADD.yml'
    shell:
        '''
        # Prescoring
        echo '## Prescored variant file' > {output.prescored};
        if [ -d {config[PrescoredFolder]} ]
        then
            for PRESCORED in $(ls {config[PrescoredFolder]}/*.tsv.gz)
            do
                cat {input} \
                | python CADD_src_scripts/extract_scored.py --header \
                    -p $PRESCORED --found_out={output.prescored}.tmp \
                > {input}.tmp;
                cat {output.prescored}.tmp >> {output.prescored}
                mv {input}.tmp {input};
            done;
            rm {output.prescored}.tmp
        fi
        mv {input} {output.novel}
        '''

rule annotation:
    input: 'CADD_tmp/{chromosome}.novel.vcf'
    output: temp('CADD_tmp/{chromosome}.anno.tsv.gz')
    conda: '../envs/CADD.yml'
    shell:
        '''
        cat {input} \
        | vep --quiet --cache --offline --dir {config[VEPpath]} \
            --buffer 1000 --no_stats --species homo_sapiens \
            --db_version={config[EnsemblDB]} --assembly {config[GenomeBuild]} \
            --format vcf --regulatory --sift b --polyphen b --per_gene --ccds --domains \
            --numbers --canonical --total_length --vcf --force_overwrite --output_file STDOUT \
        | python CADD_src_scripts/annotateVEPvcf.py \
            -c {config[ReferenceConfig]} \
        | gzip -c > {output}
        '''

rule imputation:
    input: 'CADD_tmp/{chromosome}.anno.tsv.gz'
    output: temp('CADD_tmp/{chromosome}.csv.gz')
    conda: '../envs/CADD.yml'
    shell:
        '''
        zcat {input} \
        | python CADD_src_scripts/trackTransformation.py -b \
            -c {config[ImputeConfig]} -o {output} --noheader;
        '''

rule score:
    input:
        impute='CADD_tmp/{chromosome}.csv.gz',
        anno='CADD_tmp/{chromosome}.anno.tsv.gz'
    output: temp('CADD_tmp/{chromosome}.novel.tsv')
    conda: '../envs/CADD.yml'
    shell:
        '''
        python CADD_src_scripts/predictSKmodel.py \
            -i {input.impute} -m {config[Model]} -a {input.anno} \
        | python CADD_src_scripts/max_line_hierarchy.py --all \
        | python CADD_src_scripts/appendPHREDscore.py \
            -t {config[ConversionTable]} > {output};
    
        if [ "{config[Annotation]}" = 'False' ]
        then
            cat {output} | cut -f {config[Columns]} | uniq > {output}.tmp
            mv {output}.tmp {output}
        fi
        '''

rule join:
    input:
        pre='CADD_tmp/{chromosome}.pre.tsv',
        novel='CADD_tmp/{chromosome}.novel.tsv'
    output: 'CADD/{chromosome}.tsv.gz'
    conda: '../envs/CADD.yml'
    shell:
        '''
        (
        echo "{config[Header]}";
        head -n 1 {input.novel};
        cat {input.pre} {input.novel} \
        | grep -v "^#" \
        | sort -k1,1 -k2,2n -k3,3 -k4,4 || true;
        ) | bgzip -c > {output} && tabix -s1 -b2 -e2 --force {output}
        '''

