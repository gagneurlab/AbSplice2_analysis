import pandas as pd
from tqdm import tqdm
tqdm.pandas()

def get_max(df, tps, model_name):
    df = df[df['timepoint'].isin(tps)]
    return df[model_name].max()


df = pd.read_parquet(snakemake.input['pred'])
df[['species', 'tissue_type', 'timepoint']] = df['tissue'].str.split('_', expand=True)

# get max across timepoints
tps_all = [f't{x}' for x in range(1,16)]
tps_early = [f't{x}' for x in range(1,13)]
tps_late = ['t14', 't15']

tissue = snakemake.wildcards['tissue']
model_name = snakemake.wildcards['model_name']

groupby_index = ['chrom', 'start', 'end', 'ref', 'alt', 'gene_id']

df_max_tp_all = pd.DataFrame(
    df.groupby(groupby_index).progress_apply(
        lambda df: get_max(df, tps_all, model_name))).rename(columns={0:f'{tissue}_max_all'})
df_max_tp_early = pd.DataFrame(
    df.groupby(groupby_index).progress_apply(
        lambda df: get_max(df, tps_early, model_name))).rename(columns={0:f'{tissue}_max_early'})
df_max_tp_late = pd.DataFrame(
    df.groupby(groupby_index).progress_apply(
        lambda df: get_max(df, tps_late, model_name))).rename(columns={0:f'{tissue}_max_late'})

df_max = df_max_tp_all.join(
    df_max_tp_early, how='outer').join(
    df_max_tp_late, how='outer').reset_index()

# save
df_max.to_parquet(snakemake.output['pred_max'], index=False, partition_cols=['chrom'])

