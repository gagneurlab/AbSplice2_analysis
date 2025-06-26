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
from functools import reduce

def get_abs_max_rows(df, groupby, max_col, dropna=True):
    return df.reset_index() \
        .sort_values(by=max_col, key=abs, ascending=False) \
        .drop_duplicates(subset=groupby) \
        .set_index(groupby)

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


# ======== SCRIPT STARTS HERE ========

df = pd.read_parquet(snakemake.input['pangolin_splicemap'])
df_mmsplice_splicemap = pd.read_parquet(snakemake.input['mmsplice_splicemap'])

if 'variants' in snakemake.input.keys():
    print('subsettting variants')
    df_var_subset = pd.read_parquet(snakemake.input['variants'])
    df_var_subset = df_var_subset[['chrom', 'start', 'end', 'ref', 'alt', 'gene_id']].drop_duplicates()
    df_var_subset['in_subset'] = True
    df_var_subset['variant'] = df_var_subset.apply(lambda x: f"{x['chrom']}:{x['end']}:{x['ref']}>{x['alt']}", axis=1)
    if 'variant' in df.columns:
        df = df_var_subset.set_index(['variant', 'gene_id'])[['in_subset']].join(
            df.set_index(['variant', 'gene_id']), how='left'
        ).reset_index()
    else:
        df = df_var_subset.set_index(['chrom', 'start', 'end', 'ref', 'alt', 'gene_id'])[['in_subset']].join(
            df.set_index(['chrom', 'start', 'end', 'ref', 'alt', 'gene_id']), how='left'
        ).reset_index()

    if 'variant' in df_mmsplice_splicemap.columns:
        df_mmsplice_splicemap = df_var_subset.set_index(['variant', 'gene_id'])[['in_subset']].join(
            df_mmsplice_splicemap.set_index(['variant', 'gene_id']), how='left'
        ).reset_index()
    else:
        df_mmsplice_splicemap = df_var_subset.set_index(['chrom', 'start', 'end', 'ref', 'alt', 'gene_id'])[['in_subset']].join(
            df_mmsplice_splicemap.set_index(['chrom', 'start', 'end', 'ref', 'alt', 'gene_id']), how='left'
        ).reset_index()


if 'variant' in df.columns:
    for col in ['chrom', 'start', 'end', 'ref', 'alt']:
        if col in df.columns:
            df = df.drop(col)
    df = annotate_var_compressed(df)

if df_mmsplice_splicemap.shape[0] > 0:
    if 'variant' in df_mmsplice_splicemap.columns:
        for col in ['chrom', 'start', 'end', 'ref', 'alt']:
            if col in df_mmsplice_splicemap.columns:
                df_mmsplice_splicemap = df_mmsplice_splicemap.drop(col)
        df_mmsplice_splicemap = annotate_var_compressed(df_mmsplice_splicemap)
    
    
if df_mmsplice_splicemap.empty:
    required_cols = ['chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'tissue', 'junction', 'delta_logit_psi', 'delta_psi', 'median_n', 'median_k', 'ref_psi']
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
    if 'junction' not in df_mmsplice_splicemap.columns:
        df_mmsplice_splicemap['junction'] = df_mmsplice_splicemap.apply(lambda df: f"{df['Chromosome']}:{df['junction_start']}-{df['junction_end']}:{df['junction_strand']}", axis=1)


if 'subset_cols' in snakemake.params.keys() and snakemake.params['subset_cols']:
    print('subsetting cols')
    if df_mmsplice_splicemap.shape[0] > 0:
        df_mmsplice_splicemap['j1_mmsplice'] = df_mmsplice_splicemap['junction'].apply(lambda x: x.split(':')[1].split('-')[0])
        df_mmsplice_splicemap['j2_mmsplice'] = df_mmsplice_splicemap['junction'].apply(lambda x: x.split(':')[1].split('-')[1])
        df_mmsplice_splicemap['j1_mmsplice'].astype(int)
        df_mmsplice_splicemap['j2_mmsplice'].astype(int)
    else:
        df_mmsplice_splicemap['j1_mmsplice'] = []
        df_mmsplice_splicemap['j2_mmsplice'] = []
    df_mmsplice_splicemap = df_mmsplice_splicemap[[
        'chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'tissue', 
        'ref_psi', 'median_n', 'median_k', 'delta_logit_psi', 'delta_psi', 'j1_mmsplice', 'j2_mmsplice'
    ]]
    df = df[[
        'chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'tissue', 
        'gain_score', 'gain_pos', 'loss_score', 'loss_pos',
        'pangolin_tissue_score', 'ref_psi_pangolin', 'median_n_pangolin',
    ]]
else:
    df_mmsplice_splicemap = df_mmsplice_splicemap[[
        'chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'tissue',
        'ref_psi', 'median_n', 'median_k', 'delta_logit_psi', 'delta_psi',
        'junction', 'event_type', 'splice_site'
    ]]
    df = df[[
        'chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'tissue', 
        'gain_score', 'gain_pos', 'loss_score', 'loss_pos', 
        'pangolin_tissue_score', 'ref_psi_pangolin', 'median_n_pangolin',
    ]]

    
# dropping predictions that are anyways not contributing to high scores  
if 'subset_predictions' in snakemake.params.keys() and snakemake.params['subset_predictions']:
    print('subsetting predictions')
    join_index = ['chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'tissue']
    
    if df_mmsplice_splicemap.shape[0] > 0:
        df_mmsplice_splicemap_all = df_mmsplice_splicemap.copy()
        df_mmsplice_splicemap_all['above_thresholds'] = (
            (np.abs(df_mmsplice_splicemap_all['delta_psi']) > 0.01)
            | (np.abs(df_mmsplice_splicemap_all['delta_logit_psi']) > 0.1)
        )
        df_mmsplice_splicemap = df_mmsplice_splicemap_all[df_mmsplice_splicemap_all['above_thresholds']]
        df_missing_mmsplice = df_mmsplice_splicemap_all[~df_mmsplice_splicemap_all['above_thresholds']]
        df_missing_mmsplice_max = get_abs_max_rows(df_missing_mmsplice, join_index, 'median_n').reset_index()
        df_mmsplice_splicemap = pd.concat([
            df_mmsplice_splicemap, 
            df_missing_mmsplice_max], ignore_index=True)
        df_mmsplice_splicemap = df_mmsplice_splicemap.drop(columns=['above_thresholds'])
        if 'index' in df_mmsplice_splicemap.columns:
            df_mmsplice_splicemap = df_mmsplice_splicemap.drop(columns=['index'])
        
    df_all = df.copy()
    df_all['above_thresholds'] = (
        (np.abs(df_all['loss_score']) > 0.01)
        | (np.abs(df_all['gain_score']) > 0.01)
    )
    df = df_all[df_all['above_thresholds']]
    df_missing = df_all[~df_all['above_thresholds']]
    join_index = ['chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'tissue']
    df_missing_max = get_abs_max_rows(df_missing, join_index, 'loss_score').reset_index()
    df = pd.concat([
        df, 
        df_missing_max], ignore_index=True)
    df = df.drop(columns=['above_thresholds'])
    if 'index' in df.columns:
        df = df.drop(columns=['index'])
        
# if df_mmsplice_splicemap.shape[0] > 0:
print('joining the dataframes')
join_index = ['chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'tissue']
df = df.drop_duplicates()
df_mmsplice_splicemap = df_mmsplice_splicemap.drop_duplicates()
df = df.set_index(join_index).join(
        df_mmsplice_splicemap.set_index(join_index), how='outer'
    ).reset_index()

if 'splice_site_is_expressed' not in df.columns:
    df['splice_site_is_expressed'] = (df['median_n'] > 10).astype(int)

del df_mmsplice_splicemap

def predict_in_chunks(model, df, feature_cols, chunk_size=1_000_000):
    for feature in feature_cols:
        if 'splice_site_is_expressed' in feature:
            df[feature] = df[feature].astype(int)
        else:
            df[feature] = df[feature].astype(float)
        
    preds = np.zeros(len(df), dtype=np.float32)
    print('Predicting in chunks')
    for start in range(0, len(df), chunk_size):
        end = min(start + chunk_size, len(df))
        chunk = df.iloc[start:end][feature_cols].fillna(0)
        preds[start:end] = model.predict_proba(chunk)[:, 1]

    return preds

absplice_models = snakemake.params['absplice_models']

absplice_model_names = []
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
    

def get_abs_max_rows_by_gene_chunks(df, groupby, max_col, gene_chunk_size=1000):
    all_genes = df['gene_id'].unique()
    gene_chunks = [all_genes[i:i + gene_chunk_size] for i in range(0, len(all_genes), gene_chunk_size)]
    
    chunked_results = []
    for gene_list in gene_chunks:
        chunk = df[df['gene_id'].isin(gene_list)]
        chunk = chunk.set_index(groupby_index)
        result = get_abs_max_rows(chunk, groupby, max_col)
        chunked_results.append(result)

    combined = pd.concat(chunked_results)
    final_result = get_abs_max_rows(combined, groupby, max_col)
    return final_result


if 'subset_cols' in snakemake.params.keys() and snakemake.params['subset_cols']:
    groupby_index = ['chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'tissue']

    GENE_CHUNK_SIZE = 1
    df_max = []
    all_model_names = [
        'delta_logit_psi', 'delta_psi', 'gain_score', 'loss_score',
        *absplice_model_names
    ]
    for model in tqdm([x for x in df.columns if x in all_model_names]):
        df_max.append(get_abs_max_rows_by_gene_chunks(df, groupby_index, model, gene_chunk_size=GENE_CHUNK_SIZE)[[model]].reset_index())
        # df_max.append(get_abs_max_rows_by_gene_chunks(df, groupby_index, model, gene_chunk_size=GENE_CHUNK_SIZE, n_processes=snakemake.threads)[[model]].reset_index())
    
    merged_df = reduce(
        lambda left, right: pd.merge(left, right, on=groupby_index, how='outer'),
        df_max
    )
    
    merged_df.to_parquet(
        snakemake.output['absplice_dna_pred'],
        partition_cols=['chrom'],
        index=False)
    
else:
    df.to_parquet(
        snakemake.output['absplice_dna_pred'],
        partition_cols=['chrom'])