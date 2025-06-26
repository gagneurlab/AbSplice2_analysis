import pandas as pd
import numpy as np

# read in universe
df_universe = pd.read_csv(snakemake.input['universe'])
df_universe['tissue']  = snakemake.wildcards['tissue']

# read in outliers
df_outliers = pd.read_csv(snakemake.input['outliers'])
df_outliers['outlier'] = True
df_outliers['tissue']  = snakemake.wildcards['tissue']

# join outliers to universe
df_universe = df_universe.set_index(['variant', 'gene_id', 'tissue', 'sample']).join(
    df_outliers.set_index(['variant', 'gene_id', 'tissue', 'sample']), how='left').reset_index()

# read in training data
df_preds = pd.read_parquet(snakemake.input['training_data'])

# combine
df = df_universe.set_index(['variant', 'gene_id', 'tissue']).join(
    df_preds.set_index(['variant', 'gene_id', 'tissue']), how='left').reset_index()

if snakemake.params['subset_df'] == True:
    # subset
    df = df[
        (~df['median_n_pangolin'].isna())
        | (~df['median_n_mmsplice'].isna())
    ]

    df = df[
        (np.abs(df['pangolin_score']) > 0.1)
        | (np.abs(df['delta_psi']) > 0.01)
    ]

# save
df.to_parquet(snakemake.output['training_data_combined'], index=False)