import pandas as pd
import numpy as np
from absplice.utils import get_abs_max_rows

df = pd.read_csv(snakemake.input['outliers_pval_variant_level'])

df = get_abs_max_rows(
    df.set_index(['gene_id', 'sample']), 
    ['gene_id', 'sample'], 
    'pValueGene_g_minus_log10'
).reset_index()

df.to_csv(snakemake.output['outliers_pval_gene_level'], index=False)
    