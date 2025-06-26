import pandas as pd
from absplice import SplicingOutlierResult
from absplice.utils import get_abs_max_rows

# groupby_index = ['variant', 'gene_id', 'sample', 'tissue_cat']

# df_outliers_cat = pd.concat([pd.read_csv(i) for i in snakemake.input['CAT_pval']])
# df_outliers_cat = get_abs_max_rows(
#     df_outliers_cat.set_index(groupby_index), 
#     groupby_index, 
#     'pValueGene_g_minus_log10'
# ).reset_index()

# result = SplicingOutlierResult(
#     df_absplice_dna_input = snakemake.input['absplice_dna_input'],
#     df_mmsplice_cat = snakemake.input['pred_mmsplice_cat'],
#     df_outliers_cat = df_outliers_cat,
# )

# result.absplice_rna_input.to_parquet(snakemake.output['absplice_rna_input'], partition_cols=['tissue'])

# read in files
df_absplice_dna = pd.read_parquet(snakemake.input['absplice_dna_input'])
df_mmsplice_cat = pd.read_parquet(snakemake.input['pred_mmsplice_cat']).reset_index()
df_cat_pval = pd.concat([pd.read_csv(i) for i in snakemake.input['CAT_pval']])

assert len(set(df_absplice_dna.index)) == df_absplice_dna.shape[0]

if df_mmsplice_cat.shape[0] > 0:
    # get_abs_max_rows
    df_mmsplice_cat = get_abs_max_rows(
        df_mmsplice_cat.set_index(['variant', 'gene_id', 'tissue', 'sample']), 
        groupby=['variant', 'gene_id', 'tissue', 'sample'],
        max_col='delta_psi_cat'
    )
    assert len(set(df_mmsplice_cat.index)) == df_mmsplice_cat.shape[0]
    # join
    df_absplice_rna_input = df_absplice_dna.join(
        df_mmsplice_cat[[x for x in df_mmsplice_cat.columns if 'cat' in x]], how='left'
    )
else:
    df_absplice_rna_input = df_absplice_dna.copy()

if df_cat_pval.shape[0] > 0:
    df_cat_pval = get_abs_max_rows(
        df_cat_pval.set_index(['variant', 'gene_id', 'sample']), 
        groupby=['variant', 'gene_id', 'sample'],
        max_col='pValueGene_g_minus_log10'
    )
    df_absplice_rna_input = df_absplice_rna_input.reset_index().set_index(['variant', 'gene_id', 'sample']).join(
        df_cat_pval, how='left', lsuffix='_from_cat_infer'
    )

# save
df_absplice_rna_input = df_absplice_rna_input.reset_index().set_index(['variant', 'gene_id', 'tissue', 'sample'])
df_absplice_rna_input.to_parquet(snakemake.output['absplice_rna_input'], partition_cols=['tissue'])