import pandas as pd
from functools import reduce
from absplice_scripts.utils.mapping_utils import subset_tissues


def get_abs_max_rows(df, groupby, max_col, dropna=True):
    return df.reset_index() \
        .sort_values(by=max_col, key=abs, ascending=False) \
        .drop_duplicates(subset=groupby) \
        .set_index(groupby)


absplice_models = snakemake.params['absplice_models']
groupby_index = snakemake.params['groupby_index']
groupby_index = [*groupby_index + ['Chromosome']]


df_pred = pd.read_parquet(snakemake.input['model']).reset_index()

cols_pangolin = ['variant', 'gene_id', 'tissue', 'gain_score', 'gain_pos', 'loss_score',
       'loss_pos', 'junctions_gain', 'splice_site_gain', 'ref_psi_gain',
       'median_n_gain', 'event_type_gain', 'junctions_loss',
       'splice_site_loss', 'ref_psi_loss', 'median_n_loss', 'event_type_loss',
       'pangolin_tissue_score', 'ref_psi_pangolin', 'median_n_pangolin', 'median_n_pangolin_is_na', 'splice_site_is_expressed_pangolin',
       'sample', 'Chromosome']

if 'gene_tpm' in df_pred.columns:
    cols_pangolin = [*cols_pangolin, 'gene_tpm']
if 'delta_logit_psi' in df_pred.columns:
    cols_pangolin = [*cols_pangolin, 'delta_logit_psi']
if 'delta_psi' in df_pred.columns:
    cols_pangolin = [*cols_pangolin, 'delta_psi']
if 'ref_psi' in df_pred.columns:
    cols_pangolin = [*cols_pangolin, 'ref_psi']
if 'median_n' in df_pred.columns:
    cols_pangolin = [*cols_pangolin, 'median_n']

if 'tissue_subset' in snakemake.params.keys() and snakemake.params['tissue_subset'] == True:
    df_tissue_map = pd.read_csv(snakemake.input['tissue_map'])
    tissue_map = dict(zip(df_tissue_map['tissue'], df_tissue_map['tissue_main']))
    cols_pangolin = [*cols_pangolin, 'tissue_main']
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
        
    
df_pred = df_pred[[*[x for x in cols_pangolin], *[x for x in absplice_models.keys()]]]

# Function to add suffix to all columns except those in groupby_index and model_name
def add_suffix_except_columns(df, suffix, exclude_cols):
    return df.rename(
        columns={col: f"{col}_{suffix}" if col not in exclude_cols else col for col in df.columns}
    )

print(f'GROUPBY_INDEX: {groupby_index}')
print(f'absplice_models: {absplice_models.keys()}')
# Aggregate predictions
df_pred_agg = []
for model_name, model_info in absplice_models.items():
    print('model_name:', model_name)
    # Aggregate predictions
    _df_pred_agg = get_abs_max_rows(
        df_pred.set_index(groupby_index),
        groupby=groupby_index,
        max_col=model_name
    ).reset_index()[[*cols_pangolin, model_name]]

    # Add suffix to columns (except groupby_index and model_name)
    # unique_index_cols = sorted(set([*groupby_index, 'tissue', 'tissue_main']))
    # if 'sample' in _df_pred_agg.columns:
    #     unique_index_cols = [*unique_index_cols, 'sample']
    exclude_cols = set(groupby_index + [model_name])
    _df_pred_agg = add_suffix_except_columns(_df_pred_agg, model_name, exclude_cols)

    df_pred_agg.append(_df_pred_agg)

# Merge all DataFrames in df_pred_agg on the groupby_index columns
df_pred_agg = reduce(
    lambda left, right: pd.merge(left, right, on=groupby_index, how='outer'),
    df_pred_agg
)

# print(f'unique_index_cols: {unique_index_cols}')
df_pred_agg = df_pred_agg[[
    *groupby_index, 
    *absplice_models.keys(), 
    *sorted([x for x in df_pred_agg.columns if x not in [*groupby_index, *absplice_models.keys()]])
]]

if 'tissue_subset' in snakemake.params.keys() and snakemake.params['tissue_subset'] == True:
    if 'tissue' not in df_pred_agg.columns:
        df_pred_agg['tissue'] = snakemake.wildcards['tissue_pred']

df_pred_agg.to_parquet(snakemake.output['model_agg'], partition_cols=['Chromosome'])
