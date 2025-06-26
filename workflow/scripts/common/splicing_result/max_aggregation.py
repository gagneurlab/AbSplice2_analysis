import pandas as pd
# from absplice.utils import get_abs_max_rows
from absplice_scripts.utils.mapping_utils import subset_tissues

def get_abs_max_rows(df, groupby, max_col, dropna=True):
    return df.reset_index() \
        .sort_values(by=max_col, key=abs, ascending=False) \
        .drop_duplicates(subset=groupby) \
        .set_index(groupby)

df_pred = pd.read_parquet(snakemake.input['model']).reset_index()
groupby_index = snakemake.params['groupby_index']

# if 'tissue_subset' in snakemake.params.keys() and snakemake.params['tissue_subset'] == True:
#     tissue_map = pd.read_csv(snakemake.input['tissue_map'])
#     tissue_map = dict(zip(tissue_map['tissue'], tissue_map['tissue_main']))
#     df_pred = subset_tissues(df_pred, tissue_map, chosen_tissue=snakemake.wildcards['tissue_pred'])
if 'tissue_subset' in snakemake.params.keys() and snakemake.params['tissue_subset'] == True:
    df_tissue_map = pd.read_csv(snakemake.input['tissue_map'])
    tissue_map = dict(zip(df_tissue_map['tissue'], df_tissue_map['tissue_main']))
    if snakemake.wildcards['tissue_pred'] == 'All_tissues':
        df_pred['tissue_main'] = df_pred['tissue'].map(tissue_map) # this is already done in subset_tissues
        groupby_index = [x for x in groupby_index if x != 'tissue']
    elif snakemake.wildcards['tissue_pred'] == 'All_main_tissues':
        df_pred['tissue_main'] = df_pred['tissue'].map(tissue_map) # this is already done in subset_tissues
        groupby_index[groupby_index.index('tissue')] = 'tissue_main'
    else:
        if snakemake.wildcards['tissue_pred'] in set(df_tissue_map['tissue']):
            groupby_index = groupby_index
        elif snakemake.wildcards['tissue_pred'] in set(df_tissue_map['tissue_main']):
            groupby_index[groupby_index.index('tissue')] = 'tissue_main'
        else:
            raise KeyError('%s is not in provided tissues', snakemake.wildcards['tissue_pred'])
        df_pred = subset_tissues(df_pred, tissue_map, chosen_tissue=snakemake.wildcards['tissue_pred'])


# Aggregate predictions
df_pred_agg = get_abs_max_rows(
    df_pred.set_index(groupby_index), 
    groupby=groupby_index, 
    max_col=snakemake.params['max_col']
).reset_index()

df_pred_agg['Chromosome'] = df_pred_agg['variant'].str.split(':').str[0]

df_pred_agg.to_parquet(snakemake.output['model_agg'], partition_cols=['Chromosome'])
