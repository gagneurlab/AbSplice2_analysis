import pandas as pd
import numpy as np
import pyranges as pr
import pickle
from tqdm import tqdm
from splicemap.splice_map import SpliceMap
tqdm.pandas()
import gc
import psutil

def get_pr_coords_splice_site(splice_site):
    chrom = splice_site.split(':')[0]
    pos = int(splice_site.split(':')[1])
    return chrom, pos-1, pos

def get_pr_coords_junctions(junction, position='j1'):
    chrom = junction.split(':')[0]
    if position == 'j1':
        pos = int(junction.split(':')[1].split('-')[0])
    elif position == 'j2':
        pos = int(junction.split(':')[1].split('-')[1])
    else:
        raise NotImplementedError()
    return chrom, pos-1, pos

def get_pr_coords_pangolin(row, score_type):
    chrom = row['variant'].split(':')[0]
    var_pos = int(row['variant'].split(':')[1])

    # take the variant position if the score is 0
    if np.abs(row[f'{score_type}_score']) > 0:
        pang_pos = row[f'{score_type}_pos']
        pos = var_pos + pang_pos
        return chrom, pos-1, pos
    else:
        return chrom, var_pos-1, var_pos

def pangolin_tissue_specific(row, masked=False):
    ref_psi_gain = row['ref_psi_gain']
    ref_psi_loss = row['ref_psi_loss']
    median_n_gain = row['median_n_gain']
    median_n_loss = row['median_n_loss']
    score_gain = row['gain_score']
    score_loss = row['loss_score']

    if masked:
        if ref_psi_gain > 0.8 and median_n_gain > 5:
            score_gain = 0
            ref_psi_gain = np.nan
    
        if ref_psi_loss < 0.2 or np.isnan(ref_psi_loss):
            score_loss = 0
            ref_psi_loss = np.nan
    
    if np.isnan(ref_psi_gain) and np.isnan(ref_psi_loss):
        ref_psi = np.nan
        # score = max([score_gain, score_loss], key=abs)
        score = score_gain
        median_n = median_n_gain
    elif not np.isnan(ref_psi_gain) and np.isnan(ref_psi_loss):
        ref_psi = ref_psi_gain
        score = score_gain
        median_n = median_n_gain
    elif np.isnan(ref_psi_gain) and not np.isnan(ref_psi_loss):
        ref_psi = ref_psi_loss
        score = score_loss
        median_n = median_n_loss
    else:
        if np.abs(score_gain) >= np.abs(score_loss):
            ref_psi = ref_psi_gain
            score = score_gain
            median_n = median_n_gain
        else:
            ref_psi = ref_psi_loss
            score = score_loss
            median_n = median_n_loss
            
    return score, ref_psi, median_n


def optimize_dataframe(df):
    if df.empty:
        print("DataFrame is empty. Skipping optimization.")
        return df

    print(f"Memory usage before optimization: {df.memory_usage(deep=True).sum() / 1024**3:.2f} GB")

    # Convert integers to smaller types
    for col in df.select_dtypes(include=['int64']):
        df[col] = df[col].astype('int32')

    # Convert floats to float32 and round to 4 decimal places
    for col in df.select_dtypes(include=['float64']):
        df[col] = df[col].astype('float32').round(4)

    # Convert object columns to category if unique values are small
    for col in df.select_dtypes(include=['object']):
        if len(df) > 0 and df[col].nunique() / len(df) < 0.5:
            df[col] = df[col].astype('category')

    print(f"Memory usage after optimization: {df.memory_usage(deep=True).sum() / 1024**3:.2f} GB")
    return df



SLACK = 2

# read in SpliceMaps
splicemap_cols = ['junctions', 'gene_id', 'splice_site', 'ref_psi', 'median_n']

df_splicemap = []

for splicemap_path in tqdm(snakemake.input['splicemap_5']):
    sm = SpliceMap.read_csv(splicemap_path)
    _df = sm.df
    _df = _df[splicemap_cols].copy()
    _df['tissue'] = sm.name
    _df['event_type'] = 'psi5'
    df_splicemap.append(_df)

for splicemap_path in tqdm(snakemake.input['splicemap_3']):
    sm = SpliceMap.read_csv(splicemap_path)
    _df = sm.df
    _df = _df[splicemap_cols].copy()
    _df['tissue'] = sm.name
    _df['event_type'] = 'psi3'
    df_splicemap.append(_df)

df_splicemap = pd.concat(df_splicemap)

# read in preds
df_pangolin = pd.read_parquet(snakemake.input['pred_pangolin'])
pangolin_cols = ['variant', 'gene_id', 'gain_score', 'gain_pos', 'loss_score', 'loss_pos']
df_pangolin = df_pangolin[pangolin_cols].drop_duplicates()
print(f'df_pangolin.shape: {df_pangolin.shape}')
assert df_pangolin.set_index(['variant', 'gene_id']).index.is_unique
df_pangolin = optimize_dataframe(df_pangolin)  # Add this

# Generate all unique tissue names
all_tissues = df_splicemap['tissue'].unique()
# Generate all unique variants
all_variants = df_pangolin.drop_duplicates()
# Create all variant-tissue pairs
df_variants_tissues = all_variants.merge(pd.DataFrame({'tissue': all_tissues}), how='cross')
print(f'df_variants_tissues.shape: {df_variants_tissues.shape}')
df_variants_tissues = optimize_dataframe(df_variants_tissues)  # Add this

del all_tissues
del all_variants

df_splicemap = df_splicemap[df_splicemap['gene_id'].isin(df_pangolin['gene_id'].unique())]
print(f'df_splicemap.shape: {df_splicemap.shape}')

if df_splicemap.shape[0] == 0:
    print('gene not in splicemap')
    columns = [
        'variant', 'gene_id', 'tissue', 
        'gain_score', 'gain_pos', 'loss_score', 'loss_pos']
    df = pd.DataFrame(columns=columns)
else:
    df_splicemap['Chromosome'], df_splicemap['Start'], df_splicemap['End'] = zip(*df_splicemap['splice_site'].apply(
        lambda x: get_pr_coords_splice_site(x)
    ))
    print(f'df_splicemap.shape: {df_splicemap.shape}')
    df_splicemap = optimize_dataframe(df_splicemap)  # Add this here

    df_pangolin_gain = df_pangolin.copy()
    df_pangolin_loss = df_pangolin.copy()
    del df_pangolin
    
    # get coords for the predcited splice sites
    df_pangolin_gain['Chromosome'], df_pangolin_gain['Start'], df_pangolin_gain['End'] = zip(*df_pangolin_gain.apply(
        lambda row: get_pr_coords_pangolin(row, 'gain'), axis=1
    ))
    df_pangolin_loss['Chromosome'], df_pangolin_loss['Start'], df_pangolin_loss['End'] = zip(*df_pangolin_loss.apply(
        lambda row: get_pr_coords_pangolin(row, 'loss'), axis=1
    ))
    
    # join preds with slack
    splicemap_cols = ['Chromosome', 'Start', 'End', 'junctions', 'gene_id', 'splice_site', 'ref_psi', 'median_n', 'tissue', 'event_type']
    
    df_pangolin_loss_joined = pr.PyRanges(df_pangolin_loss).join(
        pr.PyRanges(df_splicemap[splicemap_cols]), how='left', slack=SLACK
    ).df.drop(columns=['Start_b', 'End_b'])
    
    df_pangolin_gain_joined = pr.PyRanges(df_pangolin_gain).join(
        pr.PyRanges(df_splicemap[splicemap_cols]), how='left', slack=SLACK
    ).df.drop(columns=['Start_b', 'End_b'])
    
    # del unneccessary dfs
    del df_pangolin_gain
    del df_pangolin_loss
    del df_splicemap
    
    # get the correct gene ids
    df_pangolin_loss_joined = df_pangolin_loss_joined[
        (df_pangolin_loss_joined['gene_id']==df_pangolin_loss_joined['gene_id_b'])
    ].drop(columns=['gene_id_b'])
    
    df_pangolin_gain_joined = df_pangolin_gain_joined[
        (df_pangolin_gain_joined['gene_id']==df_pangolin_gain_joined['gene_id_b'])
    ].drop(columns=['gene_id_b'])
    
    if 'tissue_b' in df_pangolin_loss_joined.columns:
        df_pangolin_loss_joined = df_pangolin_loss_joined[
            (df_pangolin_loss_joined['tissue']==df_pangolin_loss_joined['tissue_b'])
        ].drop(columns=['tissue_b'])
    if 'tissue_b' in df_pangolin_gain_joined.columns:
        df_pangolin_gain_joined = df_pangolin_gain_joined[
            (df_pangolin_gain_joined['tissue']==df_pangolin_gain_joined['tissue_b'])
        ].drop(columns=['tissue_b'])
    
    # replace -1s with NaNs
    str_cols = ['junctions', 'splice_site', 'tissue', 'event_type']
    num_cols = ['ref_psi', 'median_n']
    
    df_pangolin_loss_joined[str_cols] = df_pangolin_loss_joined[str_cols].replace('-1', None)
    df_pangolin_loss_joined[num_cols] = df_pangolin_loss_joined[num_cols].replace(-1, np.nan)
    
    df_pangolin_gain_joined[str_cols] = df_pangolin_gain_joined[str_cols].replace('-1', None)
    df_pangolin_gain_joined[num_cols] = df_pangolin_gain_joined[num_cols].replace(-1, np.nan)
    
    # drop columns
    df_pangolin_gain_joined = df_pangolin_gain_joined.drop(columns=['Chromosome', 'Start', 'End'])
    df_pangolin_loss_joined = df_pangolin_loss_joined.drop(columns=['Chromosome', 'Start', 'End'])
    
    print('pangolin joining with splicemap done')
    
    # join back together loss and gain dfs
    input_index = [
        'variant',
        'gene_id',
        'tissue',
        'gain_score', 'gain_pos', 'loss_score', 'loss_pos'
    ]
    
    df = df_pangolin_gain_joined.set_index(input_index).join(
        df_pangolin_loss_joined.set_index(input_index), how='outer', lsuffix='_gain', rsuffix='_loss'
    ).reset_index()
    
    del df_pangolin_gain_joined
    del df_pangolin_loss_joined
    print(f'pangolin gain and loss joined: df.shape: {df.shape}')


# add Pangolin tissue-specific score
if df.shape[0] > 0:
    df['pangolin_tissue_score'], df['ref_psi_pangolin'], df['median_n_pangolin'] = zip(*df.apply(
        lambda row: pangolin_tissue_specific(row), axis=1
    ))
else:
    df['pangolin_tissue_score'] = []
    df['ref_psi_pangolin'] = []
    df['median_n_pangolin'] = []
    
    # Construct the new file path
    import os
    output_path_pangolin = snakemake.output['pangolin_splicemaps']
    dummy_output_path = os.path.join(
        os.path.dirname(output_path_pangolin),
        'pangolin_no_splice_sites',
        os.path.basename(output_path_pangolin).replace('parquet', 'txt')
    )

    # Make sure the parent directory exists
    os.makedirs(os.path.dirname(dummy_output_path), exist_ok=True)

    # Now you're safe to write to `output_path`
    with open(dummy_output_path, 'w') as f:
        f.write("No variant close to splice sites in provided SpliceMaps")
    

print('Before joining df to all_variants_tissues')

if 'Chromosome' in df_variants_tissues.columns:
    df_variants_tissues = df_variants_tissues.drop(columns='Chromosome')

join_index = [
    'variant', 'gene_id', 'tissue',
    'gain_score', 'gain_pos', 'loss_score', 'loss_pos']
df_variants_tissues['check'] = True

df = optimize_dataframe(df)  # Add this
df = df_variants_tissues.set_index(join_index).join(
    df.set_index(join_index)).reset_index().drop(columns='check')

df['splice_site_is_expressed_pangolin'] = (df['median_n_pangolin'] > 10).astype(int)
df['median_n_pangolin_is_na'] = (df['median_n_pangolin'].isna()).astype(int)

df = optimize_dataframe(df)  # Add this

def print_memory_usage():
    process = psutil.Process()
    print(f"Memory usage: {process.memory_info().rss / 1024**3:.2f} GB")

print_memory_usage()  # Before cleanup
gc.collect()
print_memory_usage()  # After cleanup

print(f'df.columns: {df.columns}')

if 'subset_cols' in snakemake.params.keys() and snakemake.params['subset_cols']:
    df['chrom'] = df['variant'].apply(lambda x: x.split(':')[0])
    df['pos'] = df['variant'].apply(lambda x: x.split(':')[1])
    df['pos'] = df['pos'].astype(int)
    df['ref'] = df['variant'].apply(lambda x: x.split(':')[2].split('>')[0])
    df['alt'] = df['variant'].apply(lambda x: x.split(':')[2].split('>')[1])
    df['end'] = df['pos']
    df['start'] = df['end'] - 1
    df = df[[
        'chrom', 'start', 'end', 'ref', 'alt', 'gene_id', 'tissue',
        'gain_score', 'gain_pos', 'loss_score', 'loss_pos', 'pangolin_tissue_score', 'ref_psi_pangolin', 'median_n_pangolin'
    ]]
    df = optimize_dataframe(df)
    df = df.drop_duplicates()
    df.to_parquet(
        snakemake.output['pangolin_splicemaps'],
        partition_cols=['chrom'])
    
else:
    # save
    if 'Chromosome' not in df.columns:
        df['Chromosome'] = df.apply(lambda row: row['variant'].split(':')[0], axis=1)
    df.to_parquet(
        snakemake.output['pangolin_splicemaps'],
        partition_cols=['Chromosome'])