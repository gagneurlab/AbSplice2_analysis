import pandas as pd
from tqdm import tqdm
tqdm.pandas()


df = pd.read_parquet(snakemake.input['pred'])

tissue = snakemake.wildcards['tissue']
model_name = snakemake.wildcards['model_name']

groupby_index = ['chrom', 'start', 'end', 'ref', 'alt', 'gene_id']

df_max = pd.DataFrame(
    df.groupby(groupby_index).progress_apply(
        lambda df: df[model_name].max())).rename(columns={0:f'{tissue}_max_GTEx'})

df_max = df_max.reset_index()

df_max.to_parquet(snakemake.output['pred_max'], index=False, partition_cols=['chrom'])

