import pandas as pd
from tqdm import tqdm
tqdm.pandas()

df_gtex = pd.read_parquet(snakemake.input['var_lift_gtex'])
df_dev = pd.read_parquet(snakemake.input['var_lift_dev'])

df = pd.concat([df_gtex, df_dev])
df = df.drop_duplicates()

df.to_parquet(snakemake.output['var_lift'], index=False, partition_cols='chrom')
