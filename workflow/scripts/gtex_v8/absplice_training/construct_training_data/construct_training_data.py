from splicemap.splice_map import SpliceMap
import numpy as np
import pyranges as pr
import pandas as pd
from tqdm import tqdm
tqdm.pandas()

def read_in_splicemaps(splicemap_5_list, splicemap_3_list):
    # read in SpliceMaps
    splicemap_cols = ['junctions', 'gene_id', 'splice_site', 'ref_psi', 'median_n']
    
    df_splicemap = []
    
    for splicemap_path in tqdm(splicemap_5_list):
        sm = SpliceMap.read_csv(splicemap_path)
        _df = sm.df
        _df = _df[splicemap_cols].copy()
        _df['tissue'] = sm.name
        _df['event_type'] = 'psi5'
        df_splicemap.append(_df)
    
    for splicemap_path in tqdm(splicemap_3_list):
        sm = SpliceMap.read_csv(splicemap_path)
        _df = sm.df
        _df = _df[splicemap_cols].copy()
        _df['tissue'] = sm.name
        _df['event_type'] = 'psi3'
        df_splicemap.append(_df)
    
    df_splicemap = pd.concat(df_splicemap)

    return df_splicemap


def get_pr_coords_junctions(junction, position='j1'):
    chrom = junction.split(':')[0]
    if position == 'j1':
        pos = int(junction.split(':')[1].split('-')[0])
    elif position == 'j2':
        pos = int(junction.split(':')[1].split('-')[1])
    else:
        raise NotImplementedError()
    return chrom, pos-1, pos


def get_pr_coords_pangolin(row):
    chrom = row['variant'].split(':')[0]
    var_pos = int(row['variant'].split(':')[1])

    # take the variant position if the score is 0
    if np.abs(row['pangolin_score']) > 0:
        pang_pos = row['pangolin_pos']
        pos = var_pos + pang_pos
        return chrom, pos-1, pos
    else:
        return chrom, var_pos-1, var_pos
    
    
def get_abs_max_rows(df, groupby, max_col, dropna=True):
    return df.reset_index() \
        .sort_values(by=max_col, key=abs, ascending=False) \
        .drop_duplicates(subset=groupby) \
        .set_index(groupby)
        
    
def get_abs_max_rows_mult(df, groupby, max_cols, ascending, dropna=True):
    return df.reset_index() \
        .sort_values(by=max_cols, key=abs, ascending=ascending) \
        .drop_duplicates(subset=groupby) \
        .set_index(groupby)
        
        
def get_pr_coords_splice_site(splice_site):
    chrom = splice_site.split(':')[0]
    pos = int(splice_site.split(':')[1])
    return chrom, pos-1, pos

    
def get_pangolin_tissue(df_pangolin, df_splicemap, SLACK = 2, subset_for_training=False):
    df_splicemap['Chromosome'], df_splicemap['Start'], df_splicemap['End'] = zip(*df_splicemap['splice_site'].apply(
        lambda x: get_pr_coords_splice_site(x)
    ))
    df_splicemap_all = df_splicemap.copy()

    all_tissues = df_splicemap['tissue'].unique()
    all_variants = df_pangolin.drop_duplicates()
    df_variants_tissues = all_variants.merge(pd.DataFrame({'tissue': all_tissues}), how='cross')
    
    df_pangolin_scores_loss = df_variants_tissues.drop(columns=['gain_score', 'gain_pos']).rename(columns={
        'loss_score': 'pangolin_score', 
        'loss_pos': 'pangolin_pos'
    })
    df_pangolin_scores_loss['pangolin_score_type'] = 'loss'
    df_pangolin_scores_gain = df_variants_tissues.drop(columns=['loss_score', 'loss_pos']).rename(columns={
        'gain_score': 'pangolin_score', 
        'gain_pos': 'pangolin_pos'
    })
    df_pangolin_scores_gain['pangolin_score_type'] = 'gain'
    df_pangolin_scores = pd.concat([
        df_pangolin_scores_loss,
        df_pangolin_scores_gain
    ])
    
    df_pangolin_scores['Chromosome'], df_pangolin_scores['Start'], df_pangolin_scores['End'] = zip(*df_pangolin_scores.progress_apply(
        lambda row: get_pr_coords_pangolin(row), axis=1
    ))

    df_pangolin_tissue = pr.PyRanges(df_pangolin_scores).join(
        pr.PyRanges(df_splicemap_all), how='left', slack=SLACK
    ).df.drop(columns=['Start_b', 'End_b'])

    df_pangolin_tissue = df_pangolin_tissue[
        (df_pangolin_tissue['gene_id']==df_pangolin_tissue['gene_id_b'])
        & (df_pangolin_tissue['tissue']==df_pangolin_tissue['tissue_b'])
    ].drop(columns=['gene_id_b', 'tissue_b'])

    # get only the maximum median_n for each junction
    # df_pangolin_tissue = get_abs_max_rows(df_pangolin_tissue, ['variant', 'gene_id', 'tissue', 'junctions', 'pangolin_score_type'], 'median_n').reset_index().drop(columns='index')
    if subset_for_training:
        df_pangolin_tissue = get_abs_max_rows(df_pangolin_tissue, ['variant', 'gene_id', 'tissue', 'pangolin_score_type'], 'median_n').reset_index().drop(columns='index')

    return df_pangolin_tissue, df_pangolin_scores


def annotate_median_n_categories(df_mmsplice, low_cutoff = 10, high_cutoff = 60):
    # Create the categorical column
    df_mmsplice['median_n_category'] = pd.cut(
        df_mmsplice['median_n'],
        bins=[-float('inf'), low_cutoff, high_cutoff, float('inf')],
        labels=[0, 1, 2],
        right=False
    )
    
    df_mmsplice['median_n_category'] = df_mmsplice['median_n_category'].astype(int)
    return df_mmsplice


def get_input(df_pangolin, df_mmsplice, splicemap_5_list, splicemap_3_list, subset_for_training=False):
    df_splicemap = read_in_splicemaps(splicemap_5_list, splicemap_3_list)

    df_pangolin_tissue, df_pangolin_scores = get_pangolin_tissue(df_pangolin, df_splicemap, SLACK = 2, subset_for_training=subset_for_training)
    df_pangolin_tissue = df_pangolin_tissue.drop(columns=['Chromosome', 'Start', 'End'])
    
    df_mmsplice = df_mmsplice.rename(columns={'junction': 'junctions'})
    df_mmsplice = df_mmsplice.drop(columns='gene_name')

    if subset_for_training:
        df_mmsplice = annotate_median_n_categories(df_mmsplice, low_cutoff = 10, high_cutoff = 60)
        df_mmsplice = get_abs_max_rows_mult(df_mmsplice, ['variant', 'gene_id', 'tissue'], ['median_n_category', 'delta_psi'], [False, False]).reset_index().drop(columns='index')

    # join MMSplice and Pangolin on junctions
    if subset_for_training:
        join_index = ['variant', 'gene_id', 'tissue']
    else:
        join_index = ['variant', 'gene_id', 'tissue', 'junctions']

    df_all_tissue_specific = df_pangolin_tissue.set_index(join_index).join(
        df_mmsplice.set_index(join_index), how='outer', lsuffix='_pangolin', rsuffix='_mmsplice').reset_index()

    # tissue prediction with only MMSplice info -> need to get back Pangolin part
    df_all_tissue_specific_only_mmsplice = df_all_tissue_specific[
        (df_all_tissue_specific['pangolin_score'].isna())
        & (~df_all_tissue_specific['delta_logit_psi'].isna())
    ]
    df_all_tissue_specific_only_mmsplice = df_all_tissue_specific_only_mmsplice.drop(
        columns=[x for x in df_all_tissue_specific_only_mmsplice if 'pangolin' in x])

    df_all_tissue_specific_both = df_all_tissue_specific[
        (~df_all_tissue_specific['pangolin_score'].isna())
    ]

    df_all_no_tissue_info = df_pangolin_scores.copy()
    df_all_no_tissue_info = df_all_no_tissue_info.drop(columns=['Chromosome', 'Start', 'End'])

    join_index = ['variant', 'gene_id', 'tissue']

    df_all_tissue_specific_only_mmsplice = df_all_no_tissue_info.set_index(join_index).join(
        df_all_tissue_specific_only_mmsplice.set_index(join_index), how='inner'
    ).reset_index()

    df_all_tissue_specific = pd.concat([
        df_all_tissue_specific_both,
        df_all_tissue_specific_only_mmsplice
    ])

    # __import__('pdb').set_trace()
    join_index = ['variant', 'gene_id', 'tissue', 'pangolin_score_type']
    
    tissue_indices = sorted(set(df_pangolin_scores.set_index(join_index).index).intersection(
        set(df_all_tissue_specific.set_index(join_index).index)
    ))
    
    no_tissue_info_indices = sorted(set(df_pangolin_scores.set_index(join_index).index).difference(
        set(df_all_tissue_specific.set_index(join_index).index)
    ))

    df_all = pd.concat([
        df_all_tissue_specific.set_index(join_index).loc[tissue_indices],
        df_all_no_tissue_info.set_index(join_index).loc[no_tissue_info_indices]
    ])

    return df_all


# ============ SCRIPT STARTS HERE =================

splicemap_5_list = [snakemake.input['splicemap_5']]
splicemap_3_list = [snakemake.input['splicemap_3']]

df_mmsplice = pd.read_parquet(snakemake.input['mmsplice'])
df_mmsplice = df_mmsplice[
    df_mmsplice['tissue'] == snakemake.wildcards['tissue']
]

df_pangolin = pd.read_parquet(snakemake.input['pangolin'])
df_pangolin = df_pangolin[
    df_pangolin['variant'].str.contains(f'{snakemake.wildcards["chrom"]}:')
]
df_pangolin = df_pangolin.drop(columns=['sample', 'index'])
df_pangolin = df_pangolin.drop_duplicates()

df_all_gtex_train = get_input(df_pangolin, df_mmsplice, splicemap_5_list, splicemap_3_list, subset_for_training=True)
df_all_gtex_complete = get_input(df_pangolin, df_mmsplice, splicemap_5_list, splicemap_3_list, subset_for_training=False)

df_all_gtex_train = df_all_gtex_train.reset_index()
df_all_gtex_complete = df_all_gtex_complete.reset_index()

df_all_gtex_train.to_parquet(snakemake.output['gtex_train'], index=False)
df_all_gtex_complete.to_parquet(snakemake.output['gtex_complete'], index=False)