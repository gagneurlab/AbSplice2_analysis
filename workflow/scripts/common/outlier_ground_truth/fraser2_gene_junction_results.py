import os
import pandas as pd
import numpy as np

DELTA_PSI_CUTOFF = float(snakemake.wildcards['delta_psi_cutoff'])

root_path = snakemake.input['results']

df_junc = pd.read_csv(os.path.join(root_path, 'results_per_junction.tsv'), sep='\t')
df_gene = pd.read_csv(os.path.join(root_path, 'results.tsv'), sep='\t')
# df_gene = pd.read_csv(os.path.join(root_path, 'results_gene_all.tsv'), sep='\t')

df_junc = df_junc[np.abs(df_junc['deltaPsi']) >= DELTA_PSI_CUTOFF]
df_gene = df_gene[np.abs(df_gene['deltaPsi']) >= DELTA_PSI_CUTOFF]

df_junc.to_csv(snakemake.output['junction_level'], sep='\t', index=False)
df_gene.to_csv(snakemake.output['gene_level'], sep='\t', index=False)


