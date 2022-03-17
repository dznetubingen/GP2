#!/usr/bin/env python3
import pandas as pd
import biomart 
df = pd.read_csv(snakemake.input[0],sep="\t",dtype=str)

#aggregate sample-level fields for non-comphet mode
sample_cols=['sample_id','genotype(sample,dad,mom)','depths(sample,dad,mom)', 'allele_balance(sample,dad,mom)']
others=df.loc[df['#mode']!='comphet']
sample_count=others.groupby(['chr:pos:ref:alt','#mode'])[sample_cols].agg('|'.join).reset_index()
#deal with duplicated variants (x-linked) until find a way with slivar expr
sample_count= sample_count.drop_duplicates(subset=['chr:pos:ref:alt'],keep='last')

# deal with comphet duplicated position per individaul (cuz slivar report the pair), drop not segregating comphet by counting samples in sample_id column
comphet=df.loc[df['#mode']=='comphet'].drop_duplicates()
comphet= comphet.groupby(['chr:pos:ref:alt','#mode'])[sample_cols].agg('|'.join).reset_index()
comphet_count= comphet.loc[(comphet['sample_id'].str.split("|",expand=True).count(axis=1) >= len(df['sample_id'].unique()))]

sample_df=pd.concat([sample_count,comphet_count],axis=0).sort_values(['#mode', 'chr:pos:ref:alt'])

#Extract the annotation fields (currently not good for comphet; need to find other way to merge)
#Need to update csq columns
anno=df.drop(sample_cols,axis=1).drop_duplicates(subset=['#mode','chr:pos:ref:alt'],keep='last').sort_values(['#mode', 'chr:pos:ref:alt'])

#merge back aggregate and anno fields
out=pd.merge(sample_df,anno,on=["chr:pos:ref:alt","#mode"], how='inner')
#out.drop_duplicates(['chr:pos:ref:alt'],keep='last',inplace=True)

out[['gp2_case_homalt','gp2_case_het','gp2_case_ac']] = out.GP2_affected_gt.str.split(",",expand=True)

out[['gp2_ctrl_homalt','gp2_ctrl_het','gp2_ctrl_ac']] = out.GP2_unaffected_gt.str.split(",",expand=True)

out[['fam_ctrl_homalt','fam_ctrl_het','fam_ctrl_ac']] = out.Controls.str.split(",",expand=True)

out[['fam_case_homalt','fam_case_het','fam_case_ac']] = out.Cases.str.split(",",expand=True)

def get_max_str(lst):
    return max(lst, key=len)

tmp=out.join(out.loc[:,get_max_str(out.columns)].str.split(';', expand=True).add_prefix('ann'))
tmp.drop(['GP2_affected_gt','GP2_unaffected_gt','Cases','Controls','gene',get_max_str(out.columns)], axis=1, inplace=True)

csq_column=get_max_str(out.columns).split(';')
other_column=[ c for c in tmp.columns if not(c.startswith('ann'))]
tmp.columns=other_column+csq_column

#columns to drop
allvars= tmp.drop(['gp2_case_ac','gp2_ctrl_ac','fam_ctrl_ac','fam_case_ac'],axis=1)

##need to remove amp-pd once it's fixed
allvars.rename(columns = {'Gene':'Ensembl_geneID', '#mode':'mode',
                          'NEAREST':'NEAREST_gene', 
                          'amp-pd_case_nhet-total':'amp-pd_case_nhet_3359',
                          'amp-pd_case_nhomalt-total':'amp-pd_case_nhomalt_3359',
                          'amp-pd_case_nmiss-total':'amp-pd_case_nmiss_3359'}, inplace = True)

allvars[['SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL','CADD_PHRED','gnomad_AF','gnomad_nhomalt']] = allvars[['SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL','CADD_PHRED','gnomad_AF','gnomad_nhomalt']].apply(pd.to_numeric)

## rearrange columns
column_names = ['mode','chr:pos:ref:alt','family_id','gene','Ensembl_geneID','BIOTYPE', 'transcript', 'STRAND', 'CANONICAL', 'EXON', 'Codons','Amino_acids', 'HGVSc', 'HGVSp', 'highest_impact', 'Existing_variation','ClinVar_CLNSIG','ClinVar_CLNDN','CNCR','gnomAD_pLI', 'gnomAD_oe_lof_CI90','gnomAD_oe_mis_CI90', 'gnomAD_oe_syn_CI90', 'clinvar_gene_description','MOI','CADD_RAW','CADD_PHRED','gnomad_AF', 'gnomad_popmax_af', 'gnomad_nhomalt', 'gnomad_AC','TOPMed8_AF','gp2_case_homalt', 'gp2_case_het','gp2_ctrl_homalt', 'gp2_ctrl_het', 'fam_ctrl_homalt', 'fam_ctrl_het','fam_case_homalt', 'fam_case_het','amp-pd_case_nhet_3359', 'amp-pd_case_nhomalt_3359', 'amp-pd_case_nmiss_3359','impact','NEAREST_gene', 'MAX_AF_POPS', 'SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL', 'SpliceAI_pred_SYMBOL','SpliceRegion', 'existing_InFrame_oORFs', 'existing_OutOfFrame_oORFs', 'existing_uORFs', 'five_prime_UTR_variant_annotation','five_prime_UTR_variant_consequence', 'MetaRNN_score', 'Ensembl_transcriptid', 'sample_id', 'genotype(sample,dad,mom)', 'depths(sample,dad,mom)', 'allele_balance(sample,dad,mom)','IMPACT']

allvars = allvars.reindex(columns=column_names)
allvars = allvars.sort_values(['mode','chr:pos:ref:alt','BIOTYPE','highest_impact'], ascending = (True, True, False, True))

#tier 1.1: high and moderate impact variant in protein coding genes with frq < 0.005 for dom and nhomalt <= 1
tier1= allvars[(allvars['IMPACT'].isin(['HIGH','MODERATE'])) & (allvars['gene'] != '') & ((allvars.BIOTYPE=="protein_coding") | (allvars.BIOTYPE=="retained_intron") | (allvars.BIOTYPE=="nonsense_mediated_decay"))].sort_values(['mode','chr:pos:ref:alt','BIOTYPE','highest_impact'], ascending = (True, True, False, True))

#check if there are comphets in tier1
if tier1.loc[tier1['mode']=='comphet','gene'].duplicated().any():
    c=pd.concat(g for _, g in tier1[tier1['mode']=='comphet'].groupby('gene') if len(g) > 1)
else:
    c=pd.DataFrame()

if c.shape[0] >=2:
    AD1_1 = tier1[(tier1['mode'].str.contains('dominant|denovo')) & (tier1['gnomad_AF'] < 0.005)]
    AR1_1 = tier1[(tier1['mode'].str.contains('recessive')) & (tier1['gnomad_nhomalt'] <= 1)]
    ## keep only more than 1 het in a gene for comphet
    comphet1_1=c.loc[c['gnomad_nhomalt'] <= 1]
    if  comphet1_1.empty:
        comphet1_1= pd.DataFrame()   
    elif comphet1_1.shape[0] >= 2:
        if comphet1_1['gene'].duplicated().any():
            comphet1_1=pd.concat(g for _, g in comphet1_1.groupby('gene') if len(g) > 1)
        else:
            comphet1_1=pd.DataFrame()
    else: 
        comphet1_1=pd.DataFrame()
else:
    AD1_1 = tier1[(tier1['mode'].str.contains('dominant|denovo')) & (tier1['gnomad_AF'] < 0.005)]
    AR1_1 = tier1[(tier1['mode'].str.contains('recessive')) & (tier1['gnomad_nhomalt'] <= 1)]
    comphet1_1=pd.DataFrame()

tier1_1 = pd.concat([AD1_1, AR1_1, comphet1_1], axis=0).sort_values(['mode','BIOTYPE','highest_impact','chr:pos:ref:alt'], ascending = (False, False, True,True))
tier1_2 = tier1[~tier1.apply(tuple,1).isin(tier1_1.apply(tuple,1))].sort_values(['mode','BIOTYPE','highest_impact','chr:pos:ref:alt'], ascending = (False, False, True,True))

# check comphet in tier 1_2
if tier1_2.loc[tier1_2['mode']=='comphet', 'gene'].duplicated().any():
    c=pd.concat(g for _, g in tier1_2[tier1_2['mode']=='comphet'].groupby('gene') if len(g) > 1)
else:
    c=pd.DataFrame()

if c.shape[0] >=2: 
   tier1_2=pd.concat([tier1_2[tier1_2['mode']!='comphet'],c],axis=0)
else:
   tier1_2=tier1_2[tier1_2['mode']!='comphet']

#tier 2: low impact variant in protein coding genes
tier2_tmp = allvars[((allvars['IMPACT']=='LOW') | (allvars['five_prime_UTR_variant_consequence'] !='')) & (allvars['gene'] != '') & ((allvars.BIOTYPE=="protein_coding") | (allvars.BIOTYPE=="retained_intron") | (allvars.BIOTYPE=="nonsense_mediated_decay"))].sort_values(['mode','BIOTYPE','highest_impact','chr:pos:ref:alt'], ascending = (False, False, True,True))
# add the comphet not included in tier1_1 and tier1_2
tier1_total= pd.concat([tier1_1,tier1_2],axis=0)
tier1_leftover= tier1[~tier1.apply(tuple,1).isin(tier1_total.apply(tuple,1))].sort_values(['mode','BIOTYPE','highest_impact','chr:pos:ref:alt'], ascending = (False, False, True,True))
tier2= pd.concat([tier1_leftover, tier2_tmp],axis=0).sort_values(['mode','BIOTYPE','highest_impact','chr:pos:ref:alt'], ascending = (False, False, True,True))

#check if there are comphets in tier2
if tier2.loc[tier2['mode']=='comphet', 'gene'].duplicated().any():
    c=pd.concat(g for _, g in tier2[tier2['mode']=='comphet'].groupby('gene') if len(g) > 1)
else:
    c=pd.DataFrame()

if c.shape[0] >=2:
   tier2=pd.concat([tier2[tier2['mode']!='comphet'],c],axis=0)
else:
   tier2=tier2[tier2['mode']!='comphet']

#tier3: the rest of segregating variants
s=pd.concat([tier1_1,tier1_2,tier2], axis=0)
tier3= allvars[~allvars.apply(tuple,1).isin(s.apply(tuple,1))].sort_values(['mode','BIOTYPE','highest_impact','chr:pos:ref:alt'], ascending = (False, False, True,True))

#check again the comphet
if tier2.loc[tier2['mode']=='comphet', 'gene'].duplicated().any():
    c=pd.concat(g for _, g in tier3[tier3['mode']=='comphet'].groupby('gene') if len(g) > 1)
else:
    c=pd.DataFrame()

if c.shape[0] >=2:
   tier3=pd.concat([tier3[tier3['mode']!='comphet'],c],axis=0)
else:
   tier3=tier3[tier3['mode']!='comphet']

tier1_1.to_csv(snakemake.output.impactful_1,index=False,sep="\t")
tier1_2.to_csv(snakemake.output.impactful_2,index=False,sep="\t")
tier2.to_csv(snakemake.output.low,index=False,sep="\t")
tier3.to_csv(snakemake.output.modifier,index=False,sep="\t")
