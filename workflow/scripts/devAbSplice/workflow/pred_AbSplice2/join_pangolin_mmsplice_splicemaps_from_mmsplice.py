import pandas as pd
import numpy as np
import pyranges as pr
import pickle
from tqdm import tqdm
from splicemap.splice_map import SpliceMap
tqdm.pandas()
import gc
import psutil
import os

def get_abs_max_rows(df, groupby, max_col, dropna=True):
    return df.reset_index() \
        .sort_values(by=max_col, key=abs, ascending=False) \
        .drop_duplicates(subset=groupby) \
        .set_index(groupby)

# Function to optimize DataFrame memory
def optimize_dataframe(df):
    print(f"Memory usage before optimization: {df.memory_usage(deep=True).sum() / 1024**3:.2f} GB")

    # Convert integers to smaller types
    for col in df.select_dtypes(include=['int64']):
        df[col] = df[col].astype('int32')

    # Convert floats to float32 and round to 4 decimal places
    for col in df.select_dtypes(include=['float64']):
        df[col] = df[col].astype('float32').round(4)

    # Convert object columns to category if unique values are small
    for col in df.select_dtypes(include=['object']):
        if df[col].nunique() / len(df) < 0.5:  # Only convert if cardinality is low
            df[col] = df[col].astype('category')

    print(f"Memory usage after optimization: {df.memory_usage(deep=True).sum() / 1024**3:.2f} GB")
    return df


def annotate_var_compressed(df):
    df['chrom'] = df['variant'].apply(lambda x: x.split(':')[0])
    df['pos'] = df['variant'].apply(lambda x: x.split(':')[1])
    df['pos'] = df['pos'].astype(int)
    df['ref'] = df['variant'].apply(lambda x: x.split(':')[2].split('>')[0])
    df['alt'] = df['variant'].apply(lambda x: x.split(':')[2].split('>')[1])
    df['end'] = df['pos']
    df['start'] = df['end'] - 1
    df = df.drop(columns=['variant', 'pos'])
    return df


print('reading in pangolin')
df = pd.read_parquet(snakemake.input['pangolin'])
df = optimize_dataframe(df)

print('reading in mmsplice')
df_mmsplice_splicemap = pd.read_parquet(snakemake.input['mmsplice_splicemap'])
if df_mmsplice_splicemap.empty:
    required_cols = ['variant', 'gene_id', 'tissue', 'junctions_mmsplice', 'delta_logit_psi', 'delta_psi', 'median_n', 'ref_psi']
    df_mmsplice_splicemap = pd.DataFrame(columns=required_cols)
    output_path_mmsplice = snakemake.output['absplice_dna_pred']
    dummy_output_path = os.path.join(
        os.path.dirname(output_path_mmsplice),
        'mmsplice_no_splice_sites',
        os.path.basename(output_path_mmsplice).replace('parquet', 'txt')
    )
    os.makedirs(os.path.dirname(dummy_output_path), exist_ok=True)
    with open(dummy_output_path, 'w') as f:
        f.write("No variant close to splice sites in provided SpliceMaps")   
else:
    df_mmsplice_splicemap = optimize_dataframe(df_mmsplice_splicemap)  # Add this

# add median_k from mmsplice splicemap
df_splicemap5 = []
df_splicemap3 = []
for _sm5, _sm3 in tqdm(zip(snakemake.input['splicemap_5'], snakemake.input['splicemap_3'])):
    try:
        sm5 = SpliceMap.read_csv(_sm5)
        sm3 = SpliceMap.read_csv(_sm3)
        _df5 = sm5.df
        _df3 = sm3.df
        _df5['tissue'] = sm5.name
        _df3['tissue'] = sm3.name
        df_splicemap5.append(_df5)
        df_splicemap3.append(_df3)
    except:
        pass
df_splicemap5 = pd.concat(df_splicemap5)
df_splicemap3 = pd.concat(df_splicemap3)
    
df_splicemap5 = df_splicemap5.rename(columns={'junctions': 'junction'})
df_splicemap3 = df_splicemap3.rename(columns={'junctions': 'junction'})

# read in mmsplice
df_mmsplice_splicemap5 = df_mmsplice_splicemap[df_mmsplice_splicemap['event_type'] == 'psi5']
df_mmsplice_splicemap3 = df_mmsplice_splicemap[df_mmsplice_splicemap['event_type'] == 'psi3']

# join median_k
join_index = ['junction', 'gene_id', 'tissue']

df_mmsplice_splicemap5 = df_mmsplice_splicemap5.set_index(join_index).join(
    df_splicemap5.set_index(join_index)[['median_k', 'gene_type']], how='left',
).reset_index()

df_mmsplice_splicemap3 = df_mmsplice_splicemap3.set_index(join_index).join(
    df_splicemap3.set_index(join_index)[['median_k', 'gene_type']], how='left',
).reset_index()

df_mmsplice_splicemap = pd.concat([
    df_mmsplice_splicemap5, 
    df_mmsplice_splicemap3
])

# Pangolin for all tissues
all_tissues = sorted(set(df_splicemap5['tissue']).union(set(df_splicemap3['tissue'])))
all_variants = df.drop_duplicates()
df = all_variants.merge(pd.DataFrame({'tissue': all_tissues}), how='cross')
df = optimize_dataframe(df)

if 'subset_cols' in snakemake.params.keys() and snakemake.params['subset_cols']:
    print('subsetting cols')
    df_mmsplice_splicemap['j1_mmsplice'] = df_mmsplice_splicemap['junction'].apply(lambda x: x.split(':')[1].split('-')[0])
    df_mmsplice_splicemap['j2_mmsplice'] = df_mmsplice_splicemap['junction'].apply(lambda x: x.split(':')[1].split('-')[1])
    df_mmsplice_splicemap['j1_mmsplice'].astype(int)
    df_mmsplice_splicemap['j2_mmsplice'].astype(int)
    df_mmsplice_splicemap = df_mmsplice_splicemap[[
        'variant', 'gene_id', 'tissue', 
        'ref_psi', 'median_n', 'median_k', 'delta_logit_psi', 'delta_psi', 'j1_mmsplice', 'j2_mmsplice'
    ]]
    df = df[[
        'variant', 'gene_id', 'tissue', 
        'gain_score', 'gain_pos', 'loss_score',
        'loss_pos'
    ]]
else:
    df_mmsplice_splicemap = df_mmsplice_splicemap[[
        'variant', 'gene_id', 'tissue',
        'ref_psi', 'median_n', 'median_k', 'delta_logit_psi', 'delta_psi',
        'junction', 'event_type', 'splice_site'
    ]]
    df = df[[
        'variant', 'gene_id', 'tissue', 
        'gain_score', 'gain_pos', 'loss_score', 'loss_pos', 
    ]]

    
# drop variant string column
print('compressing variant annotation')
df = annotate_var_compressed(df)
df = optimize_dataframe(df)
if df_mmsplice_splicemap.shape[0] > 0:
    df_mmsplice_splicemap = annotate_var_compressed(df_mmsplice_splicemap)
    df_mmsplice_splicemap = optimize_dataframe(df_mmsplice_splicemap)
    
# dropping predictions that are anyways not contributing to high scores  
if 'subset_predictions' in snakemake.params.keys() and snakemake.params['subset_predictions']:
    print('subsetting predictions')
    join_index = ['chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'tissue']
    
    if df_mmsplice_splicemap.shape[0] > 0:
        df_mmsplice_splicemap_all = df_mmsplice_splicemap.copy()
        df_mmsplice_splicemap = df_mmsplice_splicemap[
            (
                (np.abs(df_mmsplice_splicemap['delta_psi']) > 0.01)
                | (np.abs(df_mmsplice_splicemap['delta_logit_psi']) > 0.1)
            )
        ]
        missing_indices_mmsplice = sorted(
            set(df_mmsplice_splicemap_all.set_index(join_index).index).difference(
                set(df_mmsplice_splicemap.set_index(join_index).index)))
        df_missing_mmsplice = df_mmsplice_splicemap_all.set_index(join_index).loc[missing_indices_mmsplice].reset_index()
        df_missing_mmsplice_max = get_abs_max_rows(df_missing_mmsplice, join_index, 'delta_psi').reset_index()
        # get max for the missing indices (will anyways be small later, just for completeness)
        df_mmsplice_splicemap = pd.concat([
            df_mmsplice_splicemap, df_missing_mmsplice_max
        ])
        if 'index' in df_mmsplice_splicemap.columns:
            df_mmsplice_splicemap = df_mmsplice_splicemap.drop(columns=['index'])
        
    df_all = df.copy()
    df = df[
        (
            (np.abs(df['loss_score']) > 0.01)
            | (np.abs(df['gain_score']) > 0.01)
        )
    ]
    join_index = ['chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'tissue']
    missing_indices = sorted(
        set(df_all.set_index(join_index).index).difference(
            set(df.set_index(join_index).index)))
    df_missing = df_all.set_index(join_index).loc[missing_indices].reset_index()
    df_missing_max = get_abs_max_rows(df_missing, join_index, 'loss_score').reset_index()
    df = pd.concat([
        df, df_missing_max
    ])
    if 'index' in df.columns:
        df = df.drop(columns=['index'])
        
if df_mmsplice_splicemap.shape[0] > 0:
    print('joining the dataframes')
    join_index = ['chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'tissue']
    df = df.set_index(join_index).join(
            df_mmsplice_splicemap.set_index(join_index), how='outer'
        ).reset_index()

    print(df.columns)
    if 'splice_site_is_expressed' not in df.columns:
        df['splice_site_is_expressed'] = (df['median_n'] > 10).astype(int)

    del df_mmsplice_splicemap
    
    
if 'gene_tpm' in snakemake.input.keys():
    df_gene_tpm = pd.read_csv(snakemake.input['gene_tpm'])
    df_gene_tpm = optimize_dataframe(df_gene_tpm)  # Add this
    if 'gene_tpm' not in df.columns:
        df = df.set_index(['gene_id', 'tissue']).join(
            df_gene_tpm.set_index(['gene_id', 'tissue'])[['gene_tpm']], how='left'
        ).reset_index()

def print_memory_usage():
    process = psutil.Process()
    print(f"Memory usage: {process.memory_info().rss / 1024**3:.2f} GB")

print_memory_usage()  # Before cleanup
gc.collect()
print_memory_usage()  # After cleanup

print(f'df.columns: {df.columns}')
df = optimize_dataframe(df)

print('starting with model predictions')
# predict absplice
def predict_in_chunks(model, df, feature_cols, chunk_size=1_000_000):
    """
    Predict in memory-safe chunks.
    """
    preds = np.zeros(len(df), dtype=np.float32)
    print('Predicting in chunks')
    for start in range(0, len(df), chunk_size):
        end = min(start + chunk_size, len(df))
        chunk = df.iloc[start:end][feature_cols].fillna(0)
        preds[start:end] = model.predict_proba(chunk)[:, 1]

    return preds

absplice_models = snakemake.params['absplice_models']

absplice_model_names = []
print('right before the loop')

df['gain_score_orig_fixed'] = df['gain_score'].copy()
df['loss_score_orig_fixed'] = df['loss_score'].copy()

for model_name, model_info in absplice_models.items():
    print(f'model_name: {model_name}')
    model_path = model_info["model_path"]
    model_features = model_info["model_features"]
    with open(model_path, 'rb') as f:
        model = pickle.load(f)

    if len(set(model_features).difference(set(df.columns))) == 0:
        df[model_name] = predict_in_chunks(model, df, model_features)
        assert len(df[model_name]) == len(df)
        absplice_model_names.append(model_name)
        df = optimize_dataframe(df)
    
print('finished predictions, go to saving predictions')
print(f'df.shape: {df.shape}')

if 'subset_cols' in snakemake.params.keys() and snakemake.params['subset_cols']:
    df = df[[
        'chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'tissue',
        *absplice_model_names
    ]]
    df.to_parquet(
        snakemake.output['absplice_dna_pred'],
        partition_cols=['chrom'])
    
else:
    # save
    print('saving the predictions')
    if 'Chromosome' not in df.columns:
        df['Chromosome'] = df.apply(lambda row: row['variant'].split(':')[0], axis=1)
    df = optimize_dataframe(df)
    df.to_parquet(
        snakemake.output['absplice_dna_pred'],
        partition_cols=['Chromosome'])