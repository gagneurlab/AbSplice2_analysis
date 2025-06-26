import pandas as pd
import os
from tqdm import tqdm
tqdm.pandas()
import pickle

# read in mmsplice
df_mmsplice = pd.read_parquet(snakemake.input['mmsplice_splicemap'])
# df_mmsplice['timepoint'] = df_mmsplice['tissue'].apply(lambda x: x.split('_')[2])

# read in pangolin preds
df_pangolin = pd.read_parquet(snakemake.input['pangolin'])
df_pangolin = df_pangolin.rename(columns={
    'gain_score': 'gain_score_orig_fixed',
    'loss_score': 'loss_score_orig_fixed',
})

# join all together
join_index = ['variant', 'gene_id']
df = df_pangolin.set_index(join_index).join(df_mmsplice.set_index(join_index), how='outer').reset_index()

# predict AbSplice
model = pickle.load(open(snakemake.params['absplice_dna_model'], 'rb'))
features = 'delta_logit_psi__delta_psi__gain_score_orig_fixed__loss_score_orig_fixed__median_k__median_n'.split('__')
df_features = df[features].fillna(0).copy()
df['AbSplice_DNA'] = model.predict_proba(df_features)[:, 1]

# save results
df['Chromosome'] = df['variant'].progress_apply(lambda x: x.split(':')[0])
df.to_parquet(snakemake.output['result'], partition_cols=['Chromosome'])


