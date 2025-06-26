import pandas as pd

df = pd.concat([pd.read_parquet(i) for i in snakemake.input['training_data']])

df = df[
    (~df['median_n_pangolin'].isna())
    | (~df['median_n_mmsplice'].isna())
]

df.to_parquet(snakemake.output['training_data_subset'], index=False)