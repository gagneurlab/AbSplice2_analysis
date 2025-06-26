import pandas as pd

df = pd.read_csv(snakemake.input['SI_data'])

# subset for human
df = df[
    df['seg.id'].str.contains('hum.')
].drop(columns='Unnamed: 0')

# 1-based
df['start'] = df['start'] - 1

# pattern needs to be everything except 'n' or '-'
pattern_cols = [col for col in df.columns if col.startswith('pattern.')]
# Check for each row if all values in 'pattern.' columns are 'n' or '-'
rows_with_n_or_dash = df[pattern_cols].apply(lambda x: all(val in ['n', '-'] for val in x), axis=1)
df = df[~rows_with_n_or_dash]

df.to_csv(snakemake.output['devAS_events'], index=False)