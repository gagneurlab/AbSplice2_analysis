import pandas as pd

df = pd.read_parquet(snakemake.input['combined'])
df_lift = pd.read_parquet(snakemake.input['var_lift'])

df = df.set_index(['chrom', 'start', 'end']).join(df_lift.set_index(['chrom', 'start', 'end']), how='inner').reset_index()
df = df.rename(columns={'start': 'start_hg19', 'end': 'end_hg19'})
df['start'] = df['start_hg38'].copy()
df['end'] = df['end_hg38'].copy()

df.to_parquet(snakemake.output['combined_hg38'], index=False, partition_cols='score_category')