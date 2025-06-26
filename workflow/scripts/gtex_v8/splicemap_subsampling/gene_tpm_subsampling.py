import pandas as pd

df_gene_tpm = pd.read_csv(snakemake.input['gene_tpm'])

# Create a list to hold the new rows
new_rows = []

# Iterate through each row in the original dataframe
for index, row in tqdm(df_gene_tpm.iterrows()):
    # For each row, create 10 new rows with the modified tissue names
    for i in range(10):
        new_row = {
            'gene_id': row['gene_id'],
            'tissue': f"{row['tissue']}_4_{i}",
            'gene_tpm': row['gene_tpm']
        }
        new_rows.append(new_row)

# Convert the new rows to a dataframe and append to the original dataframe
df_new = pd.DataFrame(new_rows)

# Optionally reset the index
df_new.reset_index(drop=True, inplace=True)

df_new.to_csv(snakemake.output['gene_tpm_subsampling'], index=False)