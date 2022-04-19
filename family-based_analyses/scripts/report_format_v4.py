#!/usr/bin/env python3
import pandas as pd
import openpyxl

df = pd.read_csv(snakemake.input.tsv,sep="\t",dtype=str)
ped = pd.read_csv(snakemake.input.ped, sep='\t', header=None, names=['family_id','id','parental_id','maternal_id','sex','phenotype'])
N_ind= ped.loc[ped.phenotype==2, 'id'].count()

utr= pd.read_csv(snakemake.input.utr)
utr= utr.sort_values(['Chr','Pos'])
utr['Pos']= utr['Pos'].astype(str)

##clean utr annotation
utr['chr:pos:ref:alt']=utr[['Chr', 'Pos', 'Ref','Alt']].agg(':'.join, axis=1)
utr=utr.loc[~utr.Transcript.isnull()]
utr=utr.drop(columns=['Chr','Pos','Ref','Alt','Transcript','transcript_id'])
## get rid of all NA
utr=utr.set_index('chr:pos:ref:alt').dropna(how='all')
## keep only informative variants
cols_to_check = [c for c in utr.columns if '_gainedOrLost' in c]
cols_to_check.extend(['lost_start_codon','lost_stop_codon'])

mask = utr[cols_to_check].apply(
        lambda col:col.astype(str).str.contains(
        'gained|lost|TRUE', na=False, case=False)).any(axis=1)

utr = utr[mask].reset_index()

#aggregate sample-level fields for non-comphet mode
sample_cols=['sample_id','genotype(sample,dad,mom)','depths(sample,dad,mom)', 'allele_balance(sample,dad,mom)']
others=df.loc[df['#mode']!='comphet']
sample_count=others.groupby(['chr:pos:ref:alt','#mode'])[sample_cols].agg('|'.join).reset_index()
#deal with duplicated variants (x-linked) until find a way with slivar expr
sample_count= sample_count.drop_duplicates(subset=['chr:pos:ref:alt'],keep='last')

# deal with comphet duplicated position per individaul (cuz slivar report the pair), drop not segregating comphet by counting samples in sample_id column
comphet=df.loc[df['#mode']=='comphet'].drop_duplicates()
#drop comphet variants with highest impact as 50_intron 
comphet=comphet.loc[comphet['highest_impact']!='50_intron']
comphet= comphet.groupby(['chr:pos:ref:alt','#mode'])[sample_cols].agg('|'.join).reset_index()
comphet_count= comphet.loc[(comphet['sample_id'].str.split("|",expand=True).count(axis=1) >= len(df['sample_id'].unique()))]

sample_df=pd.concat([sample_count,comphet_count],axis=0).sort_values(['#mode', 'chr:pos:ref:alt'])

## utr output
out_utr= pd.merge(sample_df,utr, how='inner')

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

out[['amppd_case_homalt','amppd_case_het','amppd_case_ac']] = out.amp_pd_cases_gt.str.split(",",expand=True)


def get_max_str(lst):
    return max(lst, key=len)

tmp=out.join(out.loc[:,get_max_str(out.columns)].str.split(';', expand=True).add_prefix('ann'))
tmp.drop(['GP2_affected_gt','GP2_unaffected_gt','Cases','Controls','gene',get_max_str(out.columns)], axis=1, inplace=True)

csq_column=get_max_str(out.columns).split(';')
other_column=[ c for c in tmp.columns if not(c.startswith('ann'))]
tmp.columns=other_column+csq_column

##split var_synmoous
def make_OMIMlink(value):
    url = "https://www.omim.org/entry/{}"
    return '=HYPERLINK("%s", "%s")' % (url.format(value), value)

def make_Clinvarlink(value):
    url = "https://www.ncbi.nlm.nih.gov/clinvar/variation/{}"
    return '=HYPERLINK("%s", "%s")' % (url.format(value), value)


##split var_synmoous (get only clinvar accession and OMIM)
tmp['OMIM_link']=tmp["VAR_SYNONYMS"].str.extract(r'(?<=OMIM::)(.*?)(?=\.)')
tmp['OMIM_link'] = tmp['OMIM_link'].apply(lambda x: make_OMIMlink(x))
tmp['ClinVar_link']=tmp["VAR_SYNONYMS"].str.extract(r'(VCV\d*)')
tmp['ClinVar_link'] = tmp['ClinVar_link'].apply(lambda x: make_Clinvarlink(x))

tmp=tmp.astype(str)
tmp=tmp.applymap(lambda x: re.sub('.*nan.*','',x ))


#columns to drop
tmp= tmp.drop(['gp2_case_ac','gp2_ctrl_ac','fam_ctrl_ac','fam_case_ac', 'VAR_SYNONYMS'],axis=1)


##need to remove amp-pd once it's fixed
allvars.rename(columns = {'Gene':'Ensembl_geneID', '#mode':'mode',
                          'NEAREST':'NEAREST_gene', 
                          'amppd_case_homalt':'amppd_case_nhomalt_3359',
                          'amppd_case_het':'amppd_case_nhet_3359'}, inplace = True)

allvars[['SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL','CADD_PHRED','gnomad_AF','gnomad_nhomalt']] = allvars[['SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL','CADD_PHRED','gnomad_AF','gnomad_nhomalt']].apply(pd.to_numeric)

## rearrange columns
## add gene_fullname, mane, regulatroy motif, transcription factor, TSSdistance
column_names = ['mode','chr:pos:ref:alt','family_id','gene','gene_fullname','Ensembl_geneID','BIOTYPE', 'transcript','STRAND','CANONICAL','MANE_SELECT','MANE_PLUS_CLINICAL','TSL','EXON','Codons','Amino_acids','HGVSc','HGVSp','highest_impact','ClinVar_CLNSIG','ClinVar_CLNDN','Existing_variation','ClinVar_link','OMIM_link','CNCR','CADD_PHRED','gnomAD_pLI', 'gnomAD_oe_lof_CI90','gnomAD_oe_mis_CI90', 'gnomAD_oe_syn_CI90', 'clinvar_gene_description','MOI','gnomad_AF', 'gnomad_popmax_af', 'gnomad_nhomalt', 'gnomad_AC','TOPMed8_AF','gp2_case_homalt', 'gp2_case_het','gp2_ctrl_homalt', 'gp2_ctrl_het', 'fam_ctrl_homalt', 'fam_ctrl_het','fam_case_homalt', 'fam_case_het','amppd_case_nhet_3359', 'amppd_case_nhomalt_3359','impact','NEAREST_gene', 'MAX_AF_POPS','LoF','LoF_filter','LoF_flags','LoF_info','SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL', 'SpliceAI_pred_SYMBOL','SpliceRegion', 'existing_InFrame_oORFs', 'existing_OutOfFrame_oORFs', 'existing_uORFs', 'five_prime_UTR_variant_annotation','five_prime_UTR_variant_consequence', 'miRNA', 'MOTIF_NAME', 'MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE' ,'TSSDistance','TRANSCRIPTION_FACTORS','NMD','DisGeNET','MetaRNN_score', 'Ensembl_transcriptid', 'sample_id', 'genotype(sample,dad,mom)', 'depths(sample,dad,mom)', 'allele_balance(sample,dad,mom)','IMPACT']

allvars = allvars.reindex(columns=column_names)
allvars = allvars.sort_values(['mode','chr:pos:ref:alt','BIOTYPE','highest_impact'], ascending = (True, True, False, True))

##drop freq variants in gp2 due to bcftools norm 
allvars = allvars.loc[]

#tier 1.1: high and moderate impact variant in protein coding genes with frq < 0.005 for dom and nhomalt <= 1
## need to drop frequent variant in gp2 (use gp2_case_nhet-fam_case_nhet? and gp2_case_nhomalt-fam_case_nhomalt?)

tier1= allvars[(allvars['IMPACT'].isin(['HIGH','MODERATE'])) & (allvars['gene'] != '') & ((allvars.BIOTYPE=="protein_coding") | (allvars.BIOTYPE=="retained_intron") | (allvars.BIOTYPE=="nonsense_mediated_decay"))].sort_values(['mode','chr:pos:ref:alt','BIOTYPE','highest_impact'], ascending = (True, True, False, True))

AD1_1 = tier1[(tier1['mode'].str.contains('dominant|denovo')) & (tier1['gnomad_AF'] < 0.005)]
AR1_1 = tier1[(tier1['mode'].str.contains('HOM')) & (tier1['gnomad_nhomalt'] <= 1)]

#check if there are comphets in tier1
if tier1.loc[tier1['mode']=='comphet','gene'].duplicated().any():
    c=pd.concat(g for _, g in tier1[tier1['mode']=='comphet'].groupby('gene') if len(g) > 1)
else:
    c=pd.DataFrame()

##comphet in tier1 
if c.shape[0] >=2:
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
if tier3.loc[tier3['mode']=='comphet', 'gene'].duplicated().any():
    c=pd.concat(g for _, g in tier3[tier3['mode']=='comphet'].groupby('gene') if len(g) > 1)
else:
    c=pd.DataFrame()

if c.shape[0] >=2:
   tier3=pd.concat([tier3[tier3['mode']!='comphet'],c],axis=0)
else:
   tier3=tier3[tier3['mode']!='comphet']


with pd.ExcelWriter(snakemake.output.excel) as writer:
    tier1_1.to_excel(writer, sheet_name='tier_1.1',index = False, header=True)
    tier1_2.to_excel(writer, sheet_name='tier_1.2',index = False, header=True)
    tier2.to_excel(writer, sheet_name='tier_2',index = False, header=True)
    tier3.to_excel(writer, sheet_name='tier_3',index = False, header=True)
    out_utr.to_excel(writer, sheet_name='utr',index = False, header=True)
