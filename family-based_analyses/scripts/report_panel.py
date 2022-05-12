#!/usr/bin/env python3
import pandas as pd
import openpyxl
import re

df = pd.read_csv(snakemake.input.tsv,sep="\t",dtype=str)
c = pd.read_csv(snakemake.input.core_panel, header=None).iloc[:,0].to_list()
e = pd.read_csv(snakemake.input.extend_panel, header=None).iloc[:,0].to_list()
ped = pd.read_csv(snakemake.input.ped, sep='\t', header=None, names=['family_id','id','parental_id','maternal_id','sex','phenotype'])

#aggregate sample-level fields per family
sample_cols=['sample_id','genotype(sample,dad,mom)','depths(sample,dad,mom)', 'allele_balance(sample,dad,mom)']
sample_df=df.groupby(['#mode','chr:pos:ref:alt','family_id'])[sample_cols].agg('|'.join).reset_index()

#Extract the annotation fields
anno=df.drop(sample_cols,axis=1).drop_duplicates(subset=['#mode','chr:pos:ref:alt', 'family_id'],keep='last').sort_values(['#mode', 'chr:pos:ref:alt','family_id'])

#merge back aggregate and anno fields
out=pd.merge(sample_df,anno,on=['#mode','chr:pos:ref:alt','family_id'], how='inner')

def get_max_str(lst):
    return max(lst, key=len)

tmp=out.join(out.loc[:,get_max_str(out.columns)].str.split(';', expand=True).add_prefix('ann'))
csq_column=get_max_str(out.columns).split(';')
other_column=[ c for c in tmp.columns if not(c.startswith('ann'))]
tmp.columns=other_column+csq_column
tmp=tmp.loc[:,~tmp.columns.duplicated()]

def make_OMIMlink(value):
    url = "https://www.omim.org/entry/{}"
    return '=HYPERLINK("%s", "%s")' % (url.format(value), value)

def make_Clinvarlink(value):
    url = "https://www.ncbi.nlm.nih.gov/clinvar/variation/{}"
    return '=HYPERLINK("%s", "%s")' % (url.format(value), value)

##split var_synmoous (get only clinvar accession and OMIM)
tmp['OMIM_link'] = tmp['OMIM_link'].apply(lambda x: make_OMIMlink(x))
tmp['ClinVar_link']=tmp["VAR_SYNONYMS"].str.extract(r'(VCV\d*)')
tmp['ClinVar_link'] = tmp['ClinVar_link'].apply(lambda x: make_Clinvarlink(x))

tmp=tmp.astype(str)  
tmp[['OMIM_link','ClinVar_link']]=tmp[['OMIM_link','ClinVar_link']].applymap(lambda x: re.sub('.*nan.*','',x ))

#columns to drop
allvars= tmp.drop(['VAR_SYNONYMS', get_max_str(out.columns)],axis=1)

##need to remove amp-pd once it's fixed
allvars.rename(columns = {'Gene':'Ensembl_geneID', '#mode':'mode',
                          'NEAREST':'NEAREST_gene'}, inplace = True)

allvars[['SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL','CADD_PHRED','gnomad_AF','gnomad_nhomalt']] = allvars[['SpliceAI_pred_DS_AG','SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL','CADD_PHRED','gnomad_AF','gnomad_nhomalt']].apply(pd.to_numeric)

## rearrange columns
column_names = ['mode','chr:pos:ref:alt','family_id','gene','gene_fullname','Ensembl_geneID','BIOTYPE', 'transcript','STRAND','CANONICAL','MANE_SELECT','MANE_PLUS_CLINICAL','EXON','Codons','Amino_acids','HGVSc','HGVSp','highest_impact','ClinVar_CLNSIG','ClinVar_CLNDN','Existing_variation','ClinVar_link','OMIM_link','CNCR','CADD_PHRED','gnomAD_pLI', 'gnomAD_oe_lof_CI90','gnomAD_oe_mis_CI90', 'gnomAD_oe_syn_CI90', 'clinvar_gene_description','MOI','gnomad_AF', 'gnomad_popmax_af', 'gnomad_nhomalt', 'gnomad_AC','TOPMed8_AF','impact','NEAREST_gene', 'MAX_AF_POPS','LoF','LoF_filter','LoF_flags','LoF_info','SpliceAI_pred_DS_AG', 'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL', 'SpliceAI_pred_SYMBOL','SpliceRegion', 'existing_InFrame_oORFs', 'existing_OutOfFrame_oORFs', 'existing_uORFs', 'five_prime_UTR_variant_annotation','five_prime_UTR_variant_consequence', 'MetaRNN_score', 'Ensembl_transcriptid', 'sample_id', 'genotype(sample,dad,mom)', 'depths(sample,dad,mom)', 'allele_balance(sample,dad,mom)','IMPACT']

allvars = allvars.reindex(columns=column_names)
allvars = allvars.loc[allvars['BIOTYPE']=='protein_coding'].sort_values(['mode','chr:pos:ref:alt','highest_impact'], ascending = (True, True,True))

fams=list(ped['family_id'].unique())

# write out the report from each family
for family, dat in allvars.groupby('family_id'):
    with pd.ExcelWriter(f'panel/{family}.xlsx') as writer:
        dat.loc[allvars.gene.isin(c)].to_excel(writer, sheet_name='core',index = False, header=True) 
        dat.loc[allvars.gene.isin(e)].to_excel(writer, sheet_name='extend',index = False, header=True)

# write out empty report if a family doesn't have variants
fam_list=list(allvars.family_id.unique())
empty_fam=[ x for x in fams if x not in fam_list]

for family in empty_fam:
    print(family)
    c=pd.DataFrame()
    c.to_excel(f'panel/{family}.xlsx', index = False)
