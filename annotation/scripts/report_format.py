#!/usr/bin/env python3
import pandas as pd

df = pd.read_csv(snakemake.input[0],sep="\t",dtype=str)

#aggregate sample-level fields
sample_df=df.groupby(['chr:pos:ref:alt','#mode'])[['sample_id','genotype(sample,dad,mom)','depths(sample,dad,mom)', 'allele_balance(sample,dad,mom)']].agg('|'.join).reset_index()

#Extract the annotation fields (currently not good for comphet; need to find other way to merge)
anno=df.loc[df.sample_id==df.sample_id.unique()[0],['chr:pos:ref:alt','#mode', 'family_id', 'gnomad_af', 'gnomad_nhomalt', 'gnomad_ac', 'topmed_af', 'CADD_RAW', 'CADD_PHRED', 'GP2_affected_gt', 'GP2_unaffected_gt', 'amp-pd_case_nhet', 'amp-pd_case_nhomalt', 'amp-pd_case_nmiss', 'Cases', 'Controls', 'genic', 'gene', 'highest_impact', 'gene_impact_transcript_Existing_variation_MAX_AF_POPS_MetaRNN_score_Ensembl_transcriptid_SpliceAI_pred_DS_AG_SpliceAI_pred_DS_AL_SpliceAI_pred_DS_DG_SpliceAI_pred_DS_DL_SpliceAI_pred_SYMBOL_SpliceRegion_existing_InFrame_oORFs_existing_OutOfFrame_oORFs_existing_uORFs_five_prime_UTR_variant_annotation_five_prime_UTR_variant_consequence_LoFtool_Gene_BIOTYPE_STRAND_CANONICAL_NEAREST_EXON_Codons_Amino_acids_HGVSc_HGVSp_clinvar_clnsig', 'lof', 'clinvar_trait']]

#merge back aggregate and anno fields
out=pd.merge(sample_df,anno,on=["chr:pos:ref:alt","#mode"])
#out.drop_duplicates(['chr:pos:ref:alt'],keep='last',inplace=True)

out[['gp2_case_althom','gp2_case_het','gp2_case_ac']] = out.GP2_affected_gt.str.split(",",expand=True)

out[['gp2_ctrl_althom','gp2_ctrl_het','gp2_ctrl_ac']] = out.GP2_unaffected_gt.str.split(",",expand=True)

out[['fam_ctrl_althom','fam_ctrl_het','fam_ctrl_ac']] = out.Controls.str.split(",",expand=True)

out[['fam_case_althom','fam_case_het','fam_case_ac']] = out.Cases.str.split(",",expand=True)

tmp=out.join(out['gene_impact_transcript_Existing_variation_MAX_AF_POPS_MetaRNN_score_Ensembl_transcriptid_SpliceAI_pred_DS_AG_SpliceAI_pred_DS_AL_SpliceAI_pred_DS_DG_SpliceAI_pred_DS_DL_SpliceAI_pred_SYMBOL_SpliceRegion_existing_InFrame_oORFs_existing_OutOfFrame_oORFs_existing_uORFs_five_prime_UTR_variant_annotation_five_prime_UTR_variant_consequence_LoFtool_Gene_BIOTYPE_STRAND_CANONICAL_NEAREST_EXON_Codons_Amino_acids_HGVSc_HGVSp_clinvar_clnsig'].str.split('/', expand=True).add_prefix('ann'))

tmp.drop(["GP2_affected_gt","GP2_unaffected_gt","Cases","Controls","gene_impact_transcript_Existing_variation_MAX_AF_POPS_MetaRNN_score_Ensembl_transcriptid_SpliceAI_pred_DS_AG_SpliceAI_pred_DS_AL_SpliceAI_pred_DS_DG_SpliceAI_pred_DS_DL_SpliceAI_pred_SYMBOL_SpliceRegion_existing_InFrame_oORFs_existing_OutOfFrame_oORFs_existing_uORFs_five_prime_UTR_variant_annotation_five_prime_UTR_variant_consequence_LoFtool_Gene_BIOTYPE_STRAND_CANONICAL_NEAREST_EXON_Codons_Amino_acids_HGVSc_HGVSp_clinvar_clnsig"], axis=1, inplace=True)

coding=tmp.loc[(tmp.ann24!="") & (tmp.ann25!="") & (tmp.ann30.str.contains('ENST'))].reset_index(drop=True)

syn=tmp.loc[(tmp.ann24!="") & (tmp.ann26!="") & (tmp.ann29.str.contains('ENST'))].reset_index(drop=True)

utr=tmp.loc[tmp.ann28.str.contains('ENST')].reset_index(drop=True)

others=tmp.loc[tmp.ann24==""]


coding['ann32'] = coding[coding.columns[65:]].apply(\
    lambda x: ','.join(x.dropna().astype(str)),\
    axis=1)

coding=coding.iloc[:,0:66]

coding["EXON"]=coding.ann24+"|"+coding.ann25

coding["Codons"]=coding.ann26+"|"+coding.ann27

coding["Amino_acids"]=coding.ann28+"|"+coding.ann29


syn['ann31'] = syn[syn.columns[64:]].apply(\
    lambda x: ','.join(x.dropna().astype(str)),\
    axis=1)

syn=syn.iloc[:,0:66]

syn["EXON"]=syn.ann24+"|"+syn.ann25

syn["Codons"]=syn.ann26+"|"+syn.ann27

utr['ann30'] = utr[utr.columns[63:]].apply(\
    lambda x: ','.join(x.dropna().astype(str)),\
    axis=1)

utr=utr.iloc[:,0:66]

utr["EXON"]=utr.ann24+"|"+utr.ann25

##deal with clin_sig
others['ann29'] = others[others.columns[62:]].apply(lambda x: ','.join(x.dropna().astype(str)),axis=1)

others=others.iloc[:,0:66]


coding_cols=['del0','impact','transcript','Existing_variation','MAX_AF_POPS','MetaRNN_score','Ensembl_transcriptid','SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL','SpliceAI_pred_DS_DG','SpliceAI_pred_DS_DL','SpliceAI_pred_SYMBOL','SpliceRegion','existing_InFrame_oORFs','existing_OutOfFrame_oORFs','existing_uORFs','five_prime_UTR_variant_annotation','five_prime_UTR_variant_consequence','LoFtool','Gene','BIOTYPE','STRAND','CANONICAL','NEAREST','del1','del2','del3','del4','del5','del6','HGVSc','HGVSp','clinvar_clnsig','EXON','Codons','Amino_acids']

coding_col_list=list(tmp.columns[0:33])+coding_cols

coding.columns=coding_col_list

other_cols=['del0','impact','transcript','Existing_variation','MAX_AF_POPS','MetaRNN_score','Ensembl_transcriptid','SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL','SpliceAI_pred_DS_DG','SpliceAI_pred_DS_DL','SpliceAI_pred_SYMBOL','SpliceRegion','existing_InFrame_oORFs','existing_OutOfFrame_oORFs','existing_uORFs','five_prime_UTR_variant_annotation','five_prime_UTR_variant_consequence','LoFtool','Gene','BIOTYPE','STRAND','CANONICAL','NEAREST','EXON','Codons','Amino_acids','HGVSc','HGVSp','clinvar_clnsig','del4','del5','del6']

other_col_list=list(tmp.columns[0:33])+other_cols

others.columns=other_col_list

syn_cols=['del0','impact','transcript','Existing_variation','MAX_AF_POPS','MetaRNN_score','Ensembl_transcriptid','SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL','SpliceAI_pred_DS_DG','SpliceAI_pred_DS_DL','SpliceAI_pred_SYMBOL','SpliceRegion','existing_InFrame_oORFs','existing_OutOfFrame_oORFs','existing_uORFs','five_prime_UTR_variant_annotation','five_prime_UTR_variant_consequence','LoFtool','Gene','BIOTYPE','STRAND','CANONICAL','NEAREST','del1','del2','del3','del4','Amino_acids','HGVSc','HGVSp','clinvar_clnsig','del5','EXON','Codons']

syn_col_list=list(tmp.columns[0:33])+syn_cols

syn.columns=syn_col_list

utr_cols=['del0','impact','transcript','Existing_variation','MAX_AF_POPS','MetaRNN_score','Ensembl_transcriptid','SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL','SpliceAI_pred_DS_DG','SpliceAI_pred_DS_DL','SpliceAI_pred_SYMBOL','SpliceRegion','existing_InFrame_oORFs','existing_OutOfFrame_oORFs','existing_uORFs','five_prime_UTR_variant_annotation','five_prime_UTR_variant_consequence','LoFtool','Gene','BIOTYPE','STRAND','CANONICAL','NEAREST','del1','del2','Codons','Amino_acids','HGVSc','HGVSp','clinvar_clnsig','del5','del6','EXON']

utr_col_list=list(tmp.columns[0:33])+utr_cols

utr.columns=utr_col_list

coding=coding.loc[:,~coding.columns.str.startswith('del')]

utr=utr.loc[:,~utr.columns.str.startswith('del')]

syn=syn.loc[:,~syn.columns.str.startswith('del')]

others=others.loc[:,~others.columns.str.startswith('del')]

allvars=pd.concat([coding,syn,utr,others],axis=0)

allvars=allvars[['#mode','chr:pos:ref:alt','family_id','gene','Gene','BIOTYPE','STRAND','CANONICAL','transcript','EXON','Codons','Amino_acids','HGVSc','HGVSp','highest_impact','Existing_variation','clinvar_clnsig','CADD_RAW','CADD_PHRED','gnomad_af','gnomad_nhomalt','gnomad_ac','topmed_af', 'lof', 'clinvar_trait', 'gp2_case_althom', 'gp2_case_het', 'gp2_case_ac', 'amp-pd_case_nhet', 'amp-pd_case_nhomalt', 'amp-pd_case_nmiss','impact', 'NEAREST', 'MAX_AF_POPS','SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL', 'SpliceAI_pred_SYMBOL','SpliceRegion', 'existing_InFrame_oORFs', 'existing_OutOfFrame_oORFs', 'existing_uORFs', 'five_prime_UTR_variant_annotation','five_prime_UTR_variant_consequence','LoFtool','MetaRNN_score', 'Ensembl_transcriptid','sample_id','genotype(sample,dad,mom)', 'depths(sample,dad,mom)', 'allele_balance(sample,dad,mom)','genic']]

allvars[['SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL','CADD_PHRED','gnomad_af','gnomad_nhomalt']] = allvars[['SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL','CADD_PHRED','gnomad_af','gnomad_nhomalt']].apply(pd.to_numeric)

genic_tier1=allvars.loc[(allvars.genic=="genic") & (allvars.gene!="") & ((allvars.BIOTYPE=="protein_coding") | (allvars.BIOTYPE=="retained_intron") | (allvars.BIOTYPE=="nonsense_mediated_decay"))].sort_values(['#mode','BIOTYPE','highest_impact','genic'], ascending = (True, False,True,True))

genic_tier2=allvars.loc[((allvars.SpliceAI_pred_DS_AG.astype(float) >= 0.5) | (allvars.SpliceAI_pred_DS_AL.astype(float)>=0.5) | (allvars.SpliceAI_pred_DS_DG.astype(float) >=0.5)| (allvars.SpliceAI_pred_DS_DL.astype(float)>=0.5) | (allvars.CADD_PHRED >= 10)) & ~(allvars['chr:pos:ref:alt']).isin(genic_tier1['chr:pos:ref:alt'])].sort_values(['#mode','BIOTYPE','highest_impact','genic'], ascending = (True, False,True,True))

varid=list(genic_tier1['chr:pos:ref:alt'])

varid.append(list(genic_tier2['chr:pos:ref:alt']))

nongenic=allvars.loc[~allvars['chr:pos:ref:alt'].isin(varid)].sort_values(['#mode','BIOTYPE','highest_impact','genic'], ascending = (True, False,True,True))

genic_tier1.to_csv(snakemake.output.genic_tier1,index=False,sep="\t")
genic_tier2.to_csv(snakemake.output.genic_tier2,index=False,sep="\t")
nongenic.to_csv(snakemake.output.nongenic,index=False,sep="\t")
