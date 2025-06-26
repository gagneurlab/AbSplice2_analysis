import pandas as pd
from tqdm import tqdm
tqdm.pandas()

groupby_index = ['chrom', 'start', 'end', 'ref', 'alt', 'gene_id']

dev = pd.read_parquet(snakemake.input['dev']).set_index(groupby_index)
gtex = pd.read_parquet(snakemake.input['gtex']).set_index(groupby_index)

combined = dev.join(gtex, how='outer', lsuffix='_DEV', rsuffix='_GTEx').reset_index()

# get category of score
model_thresholds = {
    'max_all': {'high': 0.2, 'medium': 0.05, 'low': 0.01},
    'max_all_GTEx': {'high': 0.2, 'medium': 0.05, 'low': 0.01},
}
category_order = ['high', 'medium', 'low', 'very_low']
def get_max_category_from_models(row, model_thresholds):
    for category in category_order:
        for model, thresholds in model_thresholds.items():
            score = row.get(model, None)
            if pd.isnull(score):
                continue
            abs_score = abs(score)
            if category == 'high' and abs_score >= thresholds['high']:
                return 'high'
            elif category == 'medium' and abs_score >= thresholds['medium']:
                return 'medium'
            elif category == 'low' and abs_score >= thresholds['low']:
                return 'low'
    return 'very_low'
combined['score_category'] = combined.progress_apply(
    lambda row: get_max_category_from_models(row, model_thresholds),
    axis=1
)

# save
combined.to_parquet(snakemake.output['combined'], index=False, partition_cols=['score_category', 'chrom'])