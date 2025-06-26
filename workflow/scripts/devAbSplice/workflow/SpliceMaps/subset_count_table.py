import pandas as pd
from devAbSplice.data_kaesmann import id_map, timepoints_kaesmann

def subset_ct(ct, cols_tp):
    cols_index = ['Chromosome', 'Start', 'End', 'Strand']
    df_sm = ct[[
        *cols_index,
        *cols_tp
    ]]
    
    df_sm = df_sm.set_index(cols_index)
    
    df_sm = df_sm[~(df_sm == 0).all(axis=1)]
    return df_sm.reset_index()

species = 'human'
tissue = snakemake.wildcards['tissue']
tp = snakemake.wildcards['timepoint']
delta_t = snakemake.wildcards['delta_t']

tp_max = max([int(x.replace('t', '')) for x in timepoints_kaesmann])
tp_min = max(0, int(tp.replace('t', '')) - int(delta_t))
tp_max = min(tp_max, int(tp.replace('t', '')) + int(delta_t))

cols_tp = list()
for timepoint in range(tp_min, tp_max + 1):
    timepoint = f't{timepoint}'
    cols_tp.extend(id_map[species][tissue][timepoint])

ct = pd.read_csv(snakemake.input['raw_count_table'])
df = subset_ct(ct, cols_tp)
df.to_csv(snakemake.output['count_table'], index=False)



