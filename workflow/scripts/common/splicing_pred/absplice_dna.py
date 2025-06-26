import pandas as pd
from absplice import SplicingOutlierResult

features_dna = sorted([
    'delta_psi', 'delta_score', 'splice_site_is_expressed', 'delta_logit_psi'
])

result = SplicingOutlierResult(
    df_absplice_dna_input = snakemake.input['absplice_dna_input'],
)

result.predict_absplice_dna(
    pickle_file=snakemake.input['absplice_model'],
    median_n_cutoff=snakemake.params['median_n_cutoff'],  
    # tpm_cutoff=snakemake.params['tpm_cutoff'],
    features=features_dna, 
    abs_features=False,
    extra_info=snakemake.params['extra_info']    
)

# result._absplice_dna.to_csv(snakemake.output['absplice_dna_pred'])
# result._absplice_dna.to_parquet(snakemake.output['absplice_dna_pred'], engine='pyarrow')

df = result._absplice_dna
index = [x for x in df.index.names]
df = df.reset_index()
df['Chromosome'] = df.apply(lambda row: row['variant'].split(':')[0], axis=1)
df = df.set_index(index)
df.to_parquet(
    snakemake.output['absplice_dna_pred'],
    partition_cols=['Chromosome'])
