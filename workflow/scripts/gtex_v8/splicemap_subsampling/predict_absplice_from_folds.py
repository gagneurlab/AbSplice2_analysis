import pickle
import pandas as pd
from tqdm import tqdm
from absplice.utils import get_abs_max_rows

median_n_cutoff = 10
features = snakemake.wildcards['feature_string'].split('__')

fold = int(snakemake.wildcards['fold'])
i_range = int(snakemake.wildcards['i'])

# read in model and absplice training results
model = pickle.load(open(snakemake.input['model'], 'rb'))
df = pd.read_csv(snakemake.input['results_training'])
df_fold = df[df['fold'] == fold]
samples_fold = sorted(set(df_fold['sample']))

# read in mmsplice subsampled splicemap
df_mmsplice_splicemap_fold = []
for path in tqdm(snakemake.input['mmsplice_splicemap']):
    df_mmsplice_splicemap = pd.read_parquet(path)
    df_mmsplice_splicemap['tissue_full'] = df_mmsplice_splicemap['tissue'].copy()
    df_mmsplice_splicemap['tissue'] = df_mmsplice_splicemap['tissue_full'].apply(lambda x: '_'.join(map(str, x.split('_')[:-2])))
    df_mmsplice_splicemap['num_samples'] = df_mmsplice_splicemap['tissue_full'].apply(lambda x: int(x.split('_')[-2]))
    df_mmsplice_splicemap['i_range'] = df_mmsplice_splicemap['tissue_full'].apply(lambda x: int(x.split('_')[-1]))
    
    # subset for fold
    df_mmsplice_splicemap = df_mmsplice_splicemap[
        (df_mmsplice_splicemap['sample'].isin(samples_fold))
        & (df_mmsplice_splicemap['i_range'] == i_range)
    ]
    
    df_mmsplice_splicemap_fold.append(df_mmsplice_splicemap)
df_mmsplice_splicemap_fold = pd.concat(df_mmsplice_splicemap_fold)

# create splice_site_is_expressed feature
df_mmsplice_splicemap_fold['splice_site_is_expressed'] = (df_mmsplice_splicemap_fold['median_n'] > median_n_cutoff).astype(int)

# join mmsplice splicemap subsampled with absplice training results
join_index = ['variant', 'sample', 'tissue', 'gene_id']
df_all = df_fold.set_index(join_index).join(
    df_mmsplice_splicemap_fold.set_index(join_index), lsuffix='_from_before'
)

# perdict absplice
df_features = df_all[features].fillna(0).copy()
df_all['AbSplice_DNA_subsampled'] = model.predict_proba(df_features)[:, 1]

# max aggregate over junctions
df_all = get_abs_max_rows(df_all, groupby=join_index, max_col='AbSplice_DNA_subsampled').reset_index()

# save
df_all.to_parquet(snakemake.output['absplice'], index=False, partition_cols=['tissue'])