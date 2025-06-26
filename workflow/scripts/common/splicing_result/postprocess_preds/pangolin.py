import pandas as pd
import gzip
import pysam
import pyarrow

path = snakemake.input['model']

df = {
    'variant': [],
    'gene_id': [],
    'gain_score': [],
    'gain_pos': [],
    'loss_score': [],
    'loss_pos': [],
}

# if str(path).endswith('.gz'):
#     vcf_in =  gzip.open(path, 'rt')
# else:
#     vcf_in =  open(path, 'r')

# for line in vcf_in:
#     if not line.startswith('#') and 'Pangolin=' in line:
#         chrom = line.split('\t')[0]
#         pos = int(line.split('\t')[1])
#         ref = line.split('\t')[3]
#         alt = line.split('\t')[4]
#         per_gene_string_scores = line.split('\t')[7].split('Pangolin=')[1].split(',')
#         for gene_scores in per_gene_string_scores:
#             df['variant'].append(f'{chrom}:{pos}:{ref}>{alt}')
#             gene_id = gene_scores.split('|')[0].split('.')[0]
#             df['gene_id'].append(gene_id)
#             gain_score = float(gene_scores.split('|')[1].split(':')[1])
#             gain_pos = int(gene_scores.split('|')[1].split(':')[0])
#             loss_score = float(gene_scores.split('|')[2].split(':')[1])
#             loss_pos = int(gene_scores.split('|')[2].split(':')[0])
#             df['gain_score'].append(gain_score)
#             df['gain_pos'].append(gain_pos)
#             df['loss_score'].append(loss_score)
#             df['loss_pos'].append(loss_pos)

vcf_in = pysam.VariantFile(path, "r")

for variant in vcf_in:
    if 'Pangolin' in variant.info:
        chrom = variant.chrom
        pos = int(variant.pos)
        ref = variant.ref
        alt = variant.alts[0]
        assert variant.info['Pangolin'][0].count('ENS') == 1
        per_gene_string_scores = variant.info['Pangolin']
        for gene_scores in per_gene_string_scores:
            df['variant'].append(f'{chrom}:{pos}:{ref}>{alt}')
            gene_id = gene_scores.split('|')[0].split('.')[0]
            df['gene_id'].append(gene_id)
            gain_score = float(gene_scores.split('|')[1].split(':')[1])
            gain_pos = int(gene_scores.split('|')[1].split(':')[0])
            loss_score = float(gene_scores.split('|')[2].split(':')[1])
            loss_pos = int(gene_scores.split('|')[2].split(':')[0])
            df['gain_score'].append(gain_score)
            df['gain_pos'].append(gain_pos)
            df['loss_score'].append(loss_score)
            df['loss_pos'].append(loss_pos)

vcf_in.close()
df = pd.DataFrame(df)

# make a new column with the abs max of gain and loss scores
df['Pangolin_max_score'] = df[['gain_score', 'loss_score']].abs().max(axis=1)

if 'var_samples_df' in snakemake.input.keys():
    df_var_samples = pd.read_csv(snakemake.input['var_samples_df'])
    join_index = ['variant']
    df = df.set_index(join_index).join(
        df_var_samples.set_index(join_index), how='left'
    ).reset_index()

# df.to_parquet(snakemake.output['model_postprocess'], index=False)

df['Chromosome'] = df.apply(lambda row: row['variant'].split(':')[0], axis=1)
df.to_parquet(
    snakemake.output['model_postprocess'],
    index=False,
    partition_cols=['Chromosome'])