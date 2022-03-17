#!/usr/bin/env python
import pandas as pd

#read the input
df = pd.read_csv(snakemake.input[0],sep="\t",dtype=str)

#collapse sample-level fields
sample_df=df.groupby(['chr:pos:ref:alt','#mode'])['sample_id','genotype(sample,dad,mom)','depths(sample,dad,mom)', 'allele_balance(sample,dad,mom)'].agg('|'.join).reset_index()


#Extract the annotation fields
anno=df.loc[df.sample_id==df.sample_id.unique()[0],['chr:pos:ref:alt','#mode','gnomad_af', 'gnomad_nhomalt', 'gnomad_ac', 'topmed_af', 'CADD_RAW', 'CADD_PHRED', 'GP2_affected_gt','GP2_unaffected_gt', 'amp-pd_case_nhet', 'amp-pd_case_nhomalt','amp-pd_case_nmiss', 'Cases', 'Controls', 'gene', 'highest_impact', 'gene_impact_transcript_Existing_variation_gnomAD_genomes_controls_and_biobanks_nhomalt_MAX_AF_POPS_MetaRNN_score_SpliceAI_pred_DS_AG_SpliceAI_pred_DS_AL_SpliceAI_pred_DS_DG_SpliceAI_pred_DS_DL_BIOTYPE_STRAND_EXON_Codons_Amino_acids_HGVSc_HGVSp_clinvar_clnsig','gene_description_1', 'gene_description_2']]

out=pd.merge(sample_df,anno,on=["chr:pos:ref:alt","#mode"])
out.drop_duplicates(['chr:pos:ref:alt'],keep='last',inplace=True)
out[['gp2_case_althom','gp2_case_het','gp2_case_ac']] = out.GP2_affected_gt.str.split(",",expand=True)
out[['gp2_ctrl_althom','gp2_ctrl_het','gp2_ctrl_ac']] = out.GP2_unaffected_gt.str.split(",",expand=True)
out[['fam_ctrl_althom','fam_ctrl_het','fam_ctrl_ac']] = out.Controls.str.split(",",expand=True)
out[['fam_case_althom','fam_case_het','fam_case_ac']] = out.Cases.str.split(",",expand=True)

tmp=out.join(out['gene_impact_transcript_Existing_variation_gnomAD_genomes_controls_and_biobanks_nhomalt_MAX_AF_POPS_MetaRNN_score_SpliceAI_pred_DS_AG_SpliceAI_pred_DS_AL_SpliceAI_pred_DS_DG_SpliceAI_pred_DS_DL_BIOTYPE_STRAND_EXON_Codons_Amino_acids_HGVSc_HGVSp_clinvar_clnsig'].str.split('/', expand=True).add_prefix('ann'))
tmp.drop(["GP2_affected_gt","GP2_unaffected_gt","Cases","Controls","gene_impact_transcript_Existing_variation_gnomAD_genomes_controls_and_biobanks_nhomalt_MAX_AF_POPS_MetaRNN_score_SpliceAI_pred_DS_AG_SpliceAI_pred_DS_AL_SpliceAI_pred_DS_DG_SpliceAI_pred_DS_DL_BIOTYPE_STRAND_EXON_Codons_Amino_acids_HGVSc_HGVSp_clinvar_clnsig"], axis=1, inplace=True)

coding=tmp.loc[(tmp.ann13!="") & (tmp.ann15!="") & (tmp.ann19.str.contains('ENST'))].reset_index(drop=True)
syn=tmp.loc[(tmp.ann13!="") & (tmp.ann15!="") & (tmp.ann18.str.contains('ENST'))].reset_index(drop=True)
utr=tmp.loc[tmp.ann17.str.contains('ENST')].reset_index(drop=True)
others=tmp.loc[tmp.ann13==""]


coding["EXON"]=coding.ann13+"|"+coding.ann14
coding["Codon"]=coding.ann15+"|"+coding.ann16
coding["Amino_acids"]=coding.ann17+"|"+coding.ann18

syn["EXON"]=syn.ann13+"|"+syn.ann14
syn["Codon"]=syn.ann15+"|"+syn.ann16

utr["EXON"]=utr.ann13+"|"+utr.ann14

coding_cols=['del0','impact','transcript','Existing_variation','gnomAD_genomes_controls_and_biobanks_nhomalt','MAX_AF_POPS','MetaRNN_score','SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL','SpliceAI_pred_DS_DG','SpliceAI_pred_DS_DL','BIOTYPE','STRAND','del1','del2','del3','del4','del5','del6','HGVSc','HGVSp','clinvar_clnsig','EXON','Codons','Amino_acids']
coding_col_list=list(tmp.columns)[0:31]+coding_cols
coding.columns=coding_col_list


other_cols=['del0','impact','transcript','Existing_variation','gnomAD_genomes_controls_and_biobanks_nhomalt','MAX_AF_POPS','MetaRNN_score','SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL','SpliceAI_pred_DS_DG','SpliceAI_pred_DS_DL','BIOTYPE','STRAND','del1','del2','del3','HGVSc','HGVSp','clinvar_clnsig','del4','del5','del6']
other_col_list=list(tmp.columns)[0:31]+other_cols
others.columns=other_col_list

syn_cols=['del0','impact','transcript','Existing_variation','gnomAD_genomes_controls_and_biobanks_nhomalt','MAX_AF_POPS','MetaRNN_score','SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL','SpliceAI_pred_DS_DG','SpliceAI_pred_DS_DL','BIOTYPE','STRAND','del1','del2','del3','del4','Amino_acids','HGVSc','HGVSp','clinvar_clnsig','del5','EXON','Codons']
syn_col_list=list(tmp.columns)[0:31]+syn_cols
syn.columns=syn_col_list


utr_cols=['del0','impact','transcript','Existing_variation','gnomAD_genomes_controls_and_biobanks_nhomalt','MAX_AF_POPS','MetaRNN_score','SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL','SpliceAI_pred_DS_DG','SpliceAI_pred_DS_DL','BIOTYPE','STRAND','del1','del2','Codons','Amino_acids','HGVSc','HGVSp','clinvar_clnsig','del5','del6','EXON']
utr_col_list=list(tmp.columns)[0:31]+utr_cols
utr.columns=utr_col_list

coding=coding.loc[:,~coding.columns.str.startswith('del')]
utr=utr.loc[:,~utr.columns.str.startswith('del')]
syn=syn.loc[:,~syn.columns.str.startswith('del')]
others=others.loc[:,~others.columns.str.startswith('del')]

all=pd.concat([coding,syn,utr,others],axis=0)

all=all[['#mode','chr:pos:ref:alt','gene','BIOTYPE','highest_impact','STRAND','HGVSc', 'HGVSp','EXON', 'Codons', 'Amino_acids','gnomad_af', 'gene_description_1', 'gene_description_2','clinvar_clnsig','Existing_variation','gnomad_nhomalt', 'gnomad_ac', 'topmed_af', 'CADD_RAW', 'CADD_PHRED','amp-pd_case_nhet', 'amp-pd_case_nhomalt', 'amp-pd_case_nmiss','gp2_case_althom', 'gp2_case_het', 'gp2_case_ac', 'gp2_ctrl_althom','gp2_ctrl_het', 'gp2_ctrl_ac', 'fam_ctrl_althom', 'fam_ctrl_het','fam_ctrl_ac', 'fam_case_althom', 'fam_case_het', 'fam_case_ac','impact','transcript','gnomAD_genomes_controls_and_biobanks_nhomalt', 'MAX_AF_POPS','MetaRNN_score', 'SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL','SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL', 'sample_id', 'genotype(sample,dad,mom)', 'depths(sample,dad,mom)', 'allele_balance(sample,dad,mom)']]


all[['SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL','CADD_PHRED']] = all[['SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL','CADD_PHRED']].apply(pd.to_numeric)

protein_coding=all.loc[(all.BIOTYPE=="protein_coding") | (all.BIOTYPE=="nonsense_mediated_decay")]

code=["missense","splic","syn","stop",'frame','start','UTR']
coding_vars=protein_coding.loc[protein_coding.impact.str.contains('|'.join(code))]
splicing=protein_coding.loc[(protein_coding.SpliceAI_pred_DS_AG.astype(float) >= 0.2) | (protein_coding.SpliceAI_pred_DS_AL.astype(float)>=0.2) | (protein_coding.SpliceAI_pred_DS_DG.astype(float) >=0.2) | (protein_coding.SpliceAI_pred_DS_DL.astype(float)>=0.2) | (protein_coding.CADD_PHRED >= 10)] 


genic=pd.concat([coding_vars,splicing])\
       .sort_values('#mode')\
       .drop_duplicates(subset=['chr:pos:ref:alt'], keep='last')

genic = genic.sort_values(['#mode','BIOTYPE','highest_impact'], ascending = (True, False,True))
nongenic=all.loc[~all['chr:pos:ref:alt'].isin(genic['chr:pos:ref:alt'])]
nongenic = nongenic.sort_values(['#mode','BIOTYPE','highest_impact'], ascending = (True, False,True))


genic.to_csv("fam03_genic.tsv",index=False,sep="\t")
nongenic.to_csv("fam03_nongenic.tsv",index=False,sep="\t")
