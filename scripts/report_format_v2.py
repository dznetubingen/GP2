#!/usr/bin/env python3
import pandas as pd

df = pd.read_csv(snakemake.input[0],sep="\t",dtype=str)

#aggregate sample-level fields
sample_df=df.groupby(['chr:pos:ref:alt','#mode'])[['sample_id','genotype(sample,dad,mom)','depths(sample,dad,mom)', 'allele_balance(sample,dad,mom)']].agg('|'.join).reset_index()

#Extract the annotation fields (currently not good for comphet; need to find other way to merge)
#Need to update csq columns
anno=df.loc[df.sample_id==df.sample_id.unique()[0],['chr:pos:ref:alt','#mode', 'family_id', 'gnomad_af', 'gnomad_nhomalt', 'gnomad_ac', 'topmed_af', 'CADD_RAW', 'CADD_PHRED', 'GP2_affected_gt', 'GP2_unaffected_gt', 'amp-pd_case_nhet', 'amp-pd_case_nhomalt', 'amp-pd_case_nmiss', 'Cases', 'Controls', 'genic', 'gene', 'highest_impact', 'gene_impact_transcript;Existing_variation;MAX_AF_POPS;MetaRNN_score;Ensembl_transcriptid;SpliceAI_pred_DS_AG;SpliceAI_pred_DS_AL;SpliceAI_pred_DS_DG;SpliceAI_pred_DS_DL;SpliceAI_pred_SYMBOL;SpliceRegion;existing_InFrame_oORFs;existing_OutOfFrame_oORFs;existing_uORFs;five_prime_UTR_variant_annotation;five_prime_UTR_variant_consequence;LoFtool;Gene;BIOTYPE;STRAND;CANONICAL;NEAREST;EXON;Codons;Amino_acids;HGVSc;HGVSp;clinvar_clnsig']]

#merge back aggregate and anno fields
out=pd.merge(sample_df,anno,on=["chr:pos:ref:alt","#mode"])
#out.drop_duplicates(['chr:pos:ref:alt'],keep='last',inplace=True)

out[['gp2_case_althom','gp2_case_het','gp2_case_ac']] = out.GP2_affected_gt.str.split(",",expand=True)

out[['gp2_ctrl_althom','gp2_ctrl_het','gp2_ctrl_ac']] = out.GP2_unaffected_gt.str.split(",",expand=True)

out[['fam_ctrl_althom','fam_ctrl_het','fam_ctrl_ac']] = out.Controls.str.split(",",expand=True)

out[['fam_case_althom','fam_case_het','fam_case_ac']] = out.Cases.str.split(",",expand=True)

tmp=out.join(out['gene_impact_transcript;Existing_variation;MAX_AF_POPS;MetaRNN_score;Ensembl_transcriptid;SpliceAI_pred_DS_AG;SpliceAI_pred_DS_AL;SpliceAI_pred_DS_DG;SpliceAI_pred_DS_DL;SpliceAI_pred_SYMBOL;SpliceRegion;existing_InFrame_oORFs;existing_OutOfFrame_oORFs;existing_uORFs;five_prime_UTR_variant_annotation;five_prime_UTR_variant_consequence;LoFtool;Gene;BIOTYPE;STRAND;CANONICAL;NEAREST;EXON;Codons;Amino_acids;HGVSc;HGVSp;clinvar_clnsig'].str.split(';', expand=True).add_prefix('ann'))

tmp.drop(["GP2_affected_gt","GP2_unaffected_gt","Cases","Controls","gene_impact_transcript;Existing_variation;MAX_AF_POPS;MetaRNN_score;Ensembl_transcriptid;SpliceAI_pred_DS_AG;SpliceAI_pred_DS_AL;SpliceAI_pred_DS_DG;SpliceAI_pred_DS_DL;SpliceAI_pred_SYMBOL;SpliceRegion;existing_InFrame_oORFs;existing_OutOfFrame_oORFs;existing_uORFs;five_prime_UTR_variant_annotation;five_prime_UTR_variant_consequence;LoFtool;Gene;BIOTYPE;STRAND;CANONICAL;NEAREST;EXON;Codons;Amino_acids;HGVSc;HGVSp;clinvar_clnsig"], axis=1, inplace=True)

csq_cols=[]
cols_list=list(tmp.columns[0:33])+csq_cols
tmp.columns=cols_list

allvars=tmp[['#mode','chr:pos:ref:alt','family_id','gene','Gene','BIOTYPE','STRAND','CANONICAL','transcript','EXON','Codons','Amino_acids','HGVSc','HGVSp','highest_impact','Existing_variation','clinvar_clnsig','CADD_RAW','CADD_PHRED','gnomad_af','gnomad_nhomalt','gnomad_ac','topmed_af', 'lof', 'clinvar_trait', 'gp2_case_althom', 'gp2_case_het', 'gp2_case_ac', 'amp-pd_case_nhet', 'amp-pd_case_nhomalt', 'amp-pd_case_nmiss','impact', 'NEAREST', 'MAX_AF_POPS','SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL', 'SpliceAI_pred_SYMBOL','SpliceRegion', 'existing_InFrame_oORFs', 'existing_OutOfFrame_oORFs', 'existing_uORFs', 'five_prime_UTR_variant_annotation','five_prime_UTR_variant_consequence','LoFtool','MetaRNN_score', 'Ensembl_transcriptid','sample_id','genotype(sample,dad,mom)', 'depths(sample,dad,mom)', 'allele_balance(sample,dad,mom)','genic']]

allvars[['SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL','CADD_PHRED','gnomad_af','gnomad_nhomalt']] = allvars[['SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL','CADD_PHRED','gnomad_af','gnomad_nhomalt']].apply(pd.to_numeric)

genic_tier1=allvars.loc[(allvars.genic=="genic") & (allvars.gene!="") & ((allvars.BIOTYPE=="protein_coding") | (allvars.BIOTYPE=="retained_intron") | (allvars.BIOTYPE=="nonsense_mediated_decay"))].sort_values(['#mode','BIOTYPE','highest_impact','genic'], ascending = (True, False,True,True))

genic_tier2=allvars.loc[((allvars.SpliceAI_pred_DS_AG.astype(float) >= 0.5) | (allvars.SpliceAI_pred_DS_AL.astype(float)>=0.5) | (allvars.SpliceAI_pred_DS_DG.astype(float) >=0.5)| (allvars.SpliceAI_pred_DS_DL.astype(float)>=0.5) | (allvars.CADD_PHRED >= 10)) & ~(allvars['chr:pos:ref:alt']).isin(genic_tier1['chr:pos:ref:alt'])].sort_values(['#mode','BIOTYPE','highest_impact','genic'], ascending = (True, False,True,True))

varid=list(genic_tier1['chr:pos:ref:alt'])

varid.append(list(genic_tier2['chr:pos:ref:alt']))

nongenic=allvars.loc[~allvars['chr:pos:ref:alt'].isin(varid)].sort_values(['#mode','BIOTYPE','highest_impact','genic'], ascending = (True, False,True,True))

genic_tier1.to_csv(snakemake.output.genic_tier1,index=False,sep="\t")
genic_tier2.to_csv(snakemake.output.genic_tier2,index=False,sep="\t")
nongenic.to_csv(snakemake.output.nongenic,index=False,sep="\t")
