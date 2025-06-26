import pandas as pd
from tqdm import tqdm
import pickle
import os
import re
from absplice_scripts.visualization.benchmark import get_performance
from absplice_scripts.utils.mapping_utils import subset_tissues
from absplice.utils import get_abs_max_rows


model_dict = {
    'AbSplice_DNA_subsampled': 'AbSplice_DNA_subsampled',
    'AbSplice-DNA': 'AbSplice_DNA_absplice_dna',
}

valid_cols = [
    'outlier',
    'sample', 'gene_id', 'tissue',
    *model_dict.values(),
]

df_benchmark = list()
for i in tqdm(snakemake.input['benchmark']):
    _df = pd.read_parquet(i)
    df_benchmark.append(_df[[x for x in _df.columns if x in valid_cols]])
df_benchmark = pd.concat(df_benchmark)
assert len(set(df_benchmark.set_index(['sample', 'gene_id', 'tissue']).index)) == df_benchmark.shape[0]

# read in absplice subsampled
df_absplice_sub = list()
for i in snakemake.input['absplice_subsampled']:
    df_absplice_sub.append(pd.read_parquet(i)[['sample', 'gene_id', 'tissue', 'AbSplice_DNA_subsampled', 'fold']])
df_absplice_sub = pd.concat(df_absplice_sub)

df_absplice_sub = get_abs_max_rows(
    df_absplice_sub.set_index(['sample', 'gene_id', 'tissue']), 
    groupby=['sample', 'gene_id', 'tissue'], 
    max_col='AbSplice_DNA_subsampled')
assert len(set(df_absplice_sub.index)) == df_absplice_sub.shape[0]

# join absplice subsampled
join_index = ['sample', 'gene_id', 'tissue']
df_benchmark_all = df_benchmark.set_index(join_index).join(
    df_absplice_sub, 
    how='left', 
)
assert len(set(df_benchmark_all.index)) == df_benchmark_all.shape[0]

# performance
df_performance, performance = get_performance(
    df = df_benchmark_all, 
    outlier_column = 'outlier', 
    model_dict = model_dict,
)

# save
df_performance.to_csv(snakemake.output['df_performance'], index=False)
pickle.dump(performance, open(snakemake.output['aps_performance'], 'wb'))