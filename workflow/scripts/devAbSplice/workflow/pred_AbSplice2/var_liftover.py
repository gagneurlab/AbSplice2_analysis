import pandas as pd
from tqdm import tqdm
tqdm.pandas()

from pyliftover import LiftOver
lo = LiftOver('hg19', 'hg38')

def lift_pos(row, pos_col):
    chrom = row['chrom']
    pos = int(row[pos_col])
    try:
        pos_lifted = lo.convert_coordinate(chrom, pos)[0][1]
        return pos_lifted
    except:
        return None

df = pd.concat([pd.read_parquet(i)[['chrom', 'start', 'end']].drop_duplicates() for i in snakemake.input['absplice_dna_pred']])

df_var = df[['chrom', 'start', 'end']].drop_duplicates()

for col in ['start', 'end']:
    df_var[f'{col}_hg38'] = df_var.progress_apply(lambda x: lift_pos(x, col), axis=1)

df_var = df_var[
    (~df_var['start_hg38'].isna())
    & (~df_var['end_hg38'].isna())
]
df_var['start_hg38'] = df_var['start_hg38'].astype(int)
df_var['end_hg38'] = df_var['end_hg38'].astype(int)

df_var = df_var[
    df_var['end_hg38'] == df_var['start_hg38'] + 1
]

df_var.to_parquet(snakemake.output['var_lift'], index=False, partition_cols='chrom')

