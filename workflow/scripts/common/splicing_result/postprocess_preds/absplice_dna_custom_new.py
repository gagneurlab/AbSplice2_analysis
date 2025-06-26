import pandas as pd

df = pd.read_parquet(snakemake.input['model'])

if 'var_samples_df' in snakemake.input.keys():
    df_var_samples = pd.read_csv(snakemake.input['var_samples_df'])
    join_index = ['variant']
    df = df.set_index(join_index).join(
        df_var_samples.set_index(join_index), how='left'
    ).reset_index()
    
# save
if 'Chromosome' not in df.columns:
    df['Chromosome'] = df.apply(lambda row: row['variant'].split(':')[0], axis=1)
df.to_parquet(
    snakemake.output['model_postprocess'],
    index=False,
    partition_cols=['Chromosome'])