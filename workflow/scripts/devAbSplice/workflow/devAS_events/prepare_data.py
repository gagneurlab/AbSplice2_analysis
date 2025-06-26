import pandas as pd
import numpy as np
from tqdm import tqdm
import itertools
from multiprocessing import Pool
from devAbSplice.data_kaesmann import species_time_mapping, time_species_mapping, species_names, \
    read_species_data, get_mean_std


root_path = '/data/ouga/home/ag_gagneur/wagnern/Projects/gitlab_gagneurlab/splicing_embedding/data/raw/kaessmann/'


# read in species data
species_data = dict()
for species in tqdm(species_names):
    species_data = read_species_data(species, species_data)
    
    
# get species time data
species_data_time = dict()
for t in tqdm(time_species_mapping.keys()):
    species_data_time[t] = dict()
    for species in time_species_mapping[t].keys():
        cols = [x for x in species_data[species]['df_psi'].columns if any(pattern in x for pattern in time_species_mapping[t][species])]
        tissues = sorted(set([x.split('.')[1] for x in cols]))
        
        _df = species_data[species]['df_psi'][cols]
        if _df.shape[1] > 0:
            df_mean = get_mean_std(_df)
            species_data_time[t][species] = df_mean
        else:
            species_data_time[t][species] = pd.DataFrame()
            
            


# orthologous mapping
sup_08_orthologues = pd.read_csv(root_path + 'Supplementary_Data/Supplementary_Data_8.csv')
sup_08_orthologues = sup_08_orthologues.fillna('')
df_ortho = sup_08_orthologues.applymap(lambda x: x.split(':')[0])

# tissues = ['Brain', 'Cerebellum', 'Heart', 'Kidney', 'Liver', 'Ovary', 'Testis']

# df_ortho_time = []
# for _, row in tqdm(df_ortho.iterrows()):
    
#     for tissue in tissues:
        
#         for t in species_data_time.keys():
            
#             df_time = pd.DataFrame(row).transpose()
#             df_time['t'] = t
#             df_time['tissue'] = tissue
            
#             for species in species_data.keys():
#                 try:
#                     _df = species_data_time[t][species].loc[[row[species]]]
#                     df_time[f'{species}_mean'] = _df[f'{tissue}_mean'].values[0]
#                     df_time[f'{species}_std'] = _df[f'{tissue}_std'].values[0]
#                     df_time[f'{species}_non_null_count'] = _df[f'{tissue}_non_null_count'].values[0]
#                 except KeyError:
#                     # print(f'No orthologues for {species}')
#                     pass
#             df_ortho_time.append(df_time)

# df_ortho_time = pd.concat(df_ortho_time, axis=0, ignore_index=True)
# df_ortho_time.to_csv(snakemake.output['df_ortho_time'], index=False)

tissues = ['Brain', 'Cerebellum', 'Heart', 'Kidney', 'Liver', 'Ovary', 'Testis']

def process_row(row):
    df_ortho_time_batch = []
    index, row_data = row
    for tissue in tissues:
        for t in species_data_time.keys():
            df_time = pd.DataFrame(row_data).transpose()
            df_time['t'] = t
            df_time['tissue'] = tissue
            for species in species_data.keys():
                try:
                    _df = species_data_time[t][species].loc[[row_data[species]]]
                    df_time[f'{species}_mean'] = _df[f'{tissue}_mean'].values[0]
                    df_time[f'{species}_std'] = _df[f'{tissue}_std'].values[0]
                    df_time[f'{species}_non_null_count'] = _df[f'{tissue}_non_null_count'].values[0]
                except KeyError:
                    # print(f'No orthologues for {species}')
                    pass
            df_ortho_time_batch.append(df_time)
    return pd.concat(df_ortho_time_batch, axis=0, ignore_index=True)

def parallelize_dataframe(df, func, num_processes, batch_size):
    pool = Pool(processes=num_processes)
    results = []
    
    # Convert generator to list for batching
    df_list = list(df)
    
    # Batch processing
    num_batches = len(df_list) // batch_size
    if len(df_list) % batch_size != 0:
        num_batches += 1
    
    for i in tqdm(range(num_batches)):
        batch_start = i * batch_size
        batch_end = min((i + 1) * batch_size, len(df_list))
        batch = df_list[batch_start:batch_end]
        
        result = pool.map_async(func, batch)
        results.append(result)
    
    pool.close()
    pool.join()

    # Retrieve and concatenate results
    df_result = pd.concat([subresult for result in results for subresult in result.get()], axis=0, ignore_index=True)
    return df_result

df_ortho_time = parallelize_dataframe(df_ortho.iterrows(), process_row, num_processes=20, batch_size=1000)
df_ortho_time.to_csv(snakemake.output['df_ortho_time'], index=False)

