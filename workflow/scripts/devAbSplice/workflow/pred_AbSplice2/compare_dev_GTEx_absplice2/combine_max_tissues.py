import pandas as pd
from functools import reduce
from tqdm import tqdm
tqdm.pandas()

dfs = [pd.read_parquet(x) for x in snakemake.input['pred']]

groupby_index = ['chrom', 'start', 'end', 'ref', 'alt', 'gene_id']
joined_df = reduce(
    lambda left, right: pd.merge(left, right, on=groupby_index, how='outer'),
    dfs
)

# get max across tissues
joined_df['max_early'] = joined_df.filter(like='max_early').max(axis=1)
joined_df['max_late'] = joined_df.filter(like='max_late').max(axis=1)
joined_df['max_all'] = joined_df.filter(like='max_all').max(axis=1)

# get category of score
model_thresholds = {
    'max_all': {'high': 0.2, 'medium': 0.05, 'low': 0.01},
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
joined_df['score_category'] = joined_df.progress_apply(
    lambda row: get_max_category_from_models(row, model_thresholds),
    axis=1
)

# save
joined_df.to_parquet(snakemake.output['pred_all'], index=False, partition_cols=['score_category', 'chrom'])