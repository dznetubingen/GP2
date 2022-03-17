import pandas as pd
import glob

filenames = glob.glob( "*.xlsx")
d = []

for filename in filenames:
    xl = pd.ExcelFile(filename)
    for sheet_name in xl.sheet_names:
        df = xl.parse(sheet_name, index_col=None)
        d.append(df)
frame = pd.concat(d, axis=0, ignore_index=True)

#make multi-incident fam list (we only take variants segregating in multi-incident family))
fam=[]
flat_list = [item for sublist in fam for item in sublist]
seg=frame.loc[frame.family_id.isin(flat_list)]
hom= seg.loc[seg['mode'].str.contains('HOM')]
dom= seg.loc[seg['mode'].str.contains('dominant')]

#take all the candidate variants including singletons
hom= frame.loc[frame['mode'].str.contains('HOM')]
dom= frame.loc[frame['mode'].str.contains('dominant')]

sample_cols=['family_id','sample_id']

## recessive 
hom_count=hom.groupby(['chr:pos:ref:alt','mode'])[sample_cols].agg('|'.join).reset_index()
c=pd.concat(g for _, g in hom.groupby(['chr:pos:ref:alt','mode']) if len(g) > 1)
anno=c.groupby(['chr:pos:ref:alt','mode']).nth(0).reset_index()
not_cols=['family_id', 'sample_id', 'genotype(sample,dad,mom)', 'depths(sample,dad,mom)', 'allele_balance(sample,dad,mom)']
anno=anno[anno.columns[~anno.columns.isin(not_cols)]]
seg_hom=pd.merge(anno,hom_count,on=['chr:pos:ref:alt','mode'])
##get rid of variants with large number of carriers in AMP-PD
seg_hom[['amp-pd_case_nhet_3359','amp-pd_case_nhomalt_3359','amp-pd_case_nmiss_3359']]=seg_hom[['amp-pd_case_nhet_3359','amp-pd_case_nhomalt_3359','amp-pd_case_nmiss_3359']].apply(pd.to_numeric)
out_hom=seg_hom.loc[(seg_hom['amp-pd_case_nhomalt_3359'] < 10) & (seg_hom['highest_impact']!='67_intergenic')]


### dominant
dom_count=dom.groupby(['chr:pos:ref:alt','mode'])[sample_cols].agg('|'.join).reset_index()
c=pd.concat(g for _, g in dom.groupby(['chr:pos:ref:alt','mode']) if len(g) > 1)
anno=c.groupby(['chr:pos:ref:alt','mode']).nth(0).reset_index()
anno=anno[anno.columns[~anno.columns.isin(not_cols)]]
seg_dom=pd.merge(anno,dom_count,on=['chr:pos:ref:alt','mode'])

keep=seg_dom.loc[seg_dom['amp-pd_case_nhet_3359']=='65,9']

seg_dom[['amp-pd_case_nhet_3359','amp-pd_case_nhomalt_3359','amp-pd_case_nmiss_3359']]=seg_dom[['amp-pd_case_nhet_3359','amp-pd_case_nhomalt_3359','amp-pd_case_nmiss_3359']].apply(pd.to_numeric)
out_dom=seg_dom.loc[seg_dom['amp-pd_case_nhet_3359'] < 50 & (seg_dom['highest_impact']!='67_intergenic')]

##combine with not multi-incident familes (do not run it if not using segregating variants with multi-incident families)
others=frame.loc[~frame.family_id.isin(flat_list)]
others_count=others.groupby(['chr:pos:ref:alt','mode'])[sample_cols].agg('|'.join).reset_index()
others_count=others_count[['chr:pos:ref:alt','mode','sample_id']]

out_hom=pd.merge(seg_hom, others_count, how='inner', on=['chr:pos:ref:alt','mode'])
out_dom=pd.merge(seg_dom, others_count, how='inner', on=['chr:pos:ref:alt','mode'])

##write out
with pd.ExcelWriter('segregating_variants.xlsx') as writer:
    out_hom.to_excel(writer, sheet_name='HOM',index = False, header=True)
    out_dom.to_excel(writer, sheet_name='Dominant',index = False, header=True)
