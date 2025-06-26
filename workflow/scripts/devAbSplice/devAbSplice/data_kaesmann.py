import pandas as pd

import os
import yaml

# Get path to the config file relative to this Python file
config_path = os.path.join(os.path.dirname(__file__), '..', 'config', 'config.yaml')
config_path = os.path.abspath(config_path)

with open(config_path, 'r') as f:
    config = yaml.safe_load(f)

tissues_kaesmann = ['brain', 'cerebellum', 'heart', 'kidney', 'liver', 'ovary', 'testis']
timepoints_kaesmann = [f't{x}' for x in range(1,16)]
species_kaesmann = [
    'chicken', 
    'mouse', 
    'human',
    'macaque',
    'rat',
    'rabbit'
    #'opossum', # problem with index size, samtools index -c -m 14 {input.bamfile} {output.csi_index}
    ]

species_time_mapping = {}
species_time_mapping['human'] = {
        't1': ['.4wpc.'],
        't2': ['.5wpc.'],
        't3': ['.6wpc.'],
        't4': ['.7wpc.'],
        't5': ['.8wpc.'],
        't6': ['.9wpc.', '.10wpc.', '.11wpc.'],
        't7': ['.12wpc.'],
        't8': ['.13wpc.'],
        't9': ['.16wpc.'],
        't10': ['.18wpc.', '.19wpc.'],
        't11': ['.20wpc.'],
        't12': ['.newborn.', '.infant.', '.toddler.', ],
        't13': ['.school.', '.youngTeenager.', '.teenager.', '.oldTeenager.', '.youngAdult.'], #in youngAdult, 2 liver samples should go to t14
        't14': ['.youngMidAge.'],
        't15': ['.olderMidAge.', '.senior.', '.Senior.'],
}
species_time_mapping['macaque'] = {
        't1': [],
        't2': [],
        't3': [],
        't4': [],
        't5': [],
        't6': [],
        't7': [],
        't8': [],
        't9': [],
        't10': ['.e93.', '.e108.'],
        't11': ['.e112.', '.e123.', '.e130.'],
        't12': ['.P0.', '.P23.'],
        't13': ['.P152.', '.P183.'],
        't14': ['.P365.', '.P1095.', '.P3285.'],
        't15': ['.P5475.', '.P8030.'],
}
species_time_mapping['rat'] = {
        't1': ['.e11.'],
        't2': ['.e12.'],
        't3': ['.e13.', '.e14.'],
        't4': ['.e15.'],
        't5': ['.e16.'],
        't6': ['.e17.'],
        't7': ['.e18.'],
        't8': ['.e19.'],
        't9': ['.e20.'],
        't10': ['.P0.'],
        't11': ['.P3.'],
        't12': ['.P7.', '.P14.'],
        't13': ['.P42.'],
        't14': ['.P112.'],
        't15': [],
}
species_time_mapping['mouse'] = {
        't1': ['.e10.5.'],
        't2': ['.e11.5.'],
        't3': ['.e12.5.'],
        't4': ['.e13.5.'],
        't5': ['.e14.5.'],
        't6': ['.e15.5.'],
        't7': ['.e16.5.'],
        't8': ['.e17.5.'],
        't9': ['.e18.5.'],
        't10': ['.P0.'],
        't11': ['.P3.'],
        't12': ['.P14.'],
        't13': ['.P28.'],
        't14': ['.P63.'],
        't15': [],
}
species_time_mapping['rabbit'] = {
        't1': ['.e12.'],
        't2': ['.e13.'],
        't3': ['.e14.', '.e15.5.'],
        't4': ['.e16.5.'],
        't5': ['.e18.'],
        't6': ['.e19.5.'],
        't7': ['.e21.'],
        't8': ['.e23.'],
        't9': ['.e24.'],
        't10': ['.e27.'],
        't11': ['.P0.'],
        't12': ['.P14.'],
        't13': ['.P84.'],
        't14': ['.P186.'], #P180 not there, P186 contains info from t15 as well
        't15': [],
}
species_time_mapping['opossum'] = {
        't1': ['.13.5.'], # corresponds to e13.5
        't2': ['.14.'], # corresponds to P0
        't3': ['.16.'], # corresponds to P2
        't4': ['.16.'], # corresponds to P2
        't5': ['.18.'], # corresponds to P4
        't6': ['.18.'], # corresponds to P4
        't7': ['.20.'], # corresponds to P6
        't8': ['.24.'], # corresponds to P10
        't9': ['.28.'], # corresponds to P14
        't10': ['.35.'], # corresponds to P21
        't11': ['.42.'], # corresponds to P28
        't12': ['.56.', '.74.'], # corresponds to P42, P60
        't13': ['.104.'], # corresponds to P90
        't14': ['.134.'], # corresponds to P120
        't15': ['.164.', '.194.'], # corresponds to P270, P180, P420
}
species_time_mapping['chicken'] = {
        't1': [],
        't2': [],
        't3': [],
        't4': [],
        't5': [],
        't6': [],
        't7': ['.e10.'],
        't8': ['.e12.'],
        't9': ['.e12.'],
        't10': ['.e12.'],
        't11': ['.e14.'],
        't12': ['.e17.', '.P0.'],
        't13': ['.P7.', '.P35.', '.P70.'],
        't14': ['.P155.'],
        't15': []
}
species_names = species_time_mapping.keys()


time_species_mapping = {}

for species, time_data in species_time_mapping.items():
    for timepoint, value in time_data.items():
        if timepoint in time_species_mapping:
            time_species_mapping[timepoint][species] = value
        else:
            time_species_mapping[timepoint] = {species: value}
            
            
DATA_DIR = config['DATA_DIR']

column_names_sajr = [
    'chrom',
    'annotation_source',
    'type',
    'start',
    'end',
    'not_known',
    'strand',
    'not_known2',
    'annotation'  
]
def read_species_data(species_name, species_data):
    species_data[species_name] = dict()
    species_data[species_name]['df_e'] = pd.read_csv(DATA_DIR + f'{species_name}.e.gz').rename(columns={'Unnamed: 0': 'annotation'}).set_index('annotation')
    species_data[species_name]['df_i'] = pd.read_csv(DATA_DIR + f'{species_name}.i.gz').rename(columns={'Unnamed: 0': 'annotation'}).set_index('annotation')
    species_data[species_name]['df_dev'] = pd.read_csv(DATA_DIR + f'{species_name}.devAS.gz').rename(columns={'Unnamed: 0': 'annotation'}).set_index('annotation')
    species_data[species_name]['df_psi'] = pd.read_csv(DATA_DIR + f'{species_name}.psi.gz').rename(columns={'Unnamed: 0': 'annotation'}).set_index('annotation')
    species_data[species_name]['df_sajr'] = pd.read_csv(DATA_DIR + f'{species_name}.sajr.gz', sep='\t', skiprows=4, names=column_names_sajr)
    return species_data

tissues = ['Brain', 'Cerebellum', 'Heart', 'Kidney', 'Liver', 'Ovary', 'Testis']
def get_mean_std(_df):
    df_mean = list()
    for tissue in tissues:
        col_tissue = [x for x in _df.columns if x.split('.')[1] == tissue]
        mean_values = _df[col_tissue].mean(skipna=True, axis=1)
        std_values = _df[col_tissue].std(skipna=True, ddof=1, axis=1)
        non_null_count = _df[col_tissue].count(axis=1)
        df_mean.append(pd.DataFrame({
            f'{tissue}_mean': mean_values,
            f'{tissue}_std': std_values,
            f'{tissue}_non_null_count': non_null_count
        }))
    df_mean = pd.concat(df_mean, axis=1)
    return df_mean



# Species	Genome assembly	Ensembl annotation versions
# Gallus gallus	Galgal4	v84
# Homo sapiens	GRCh37.73	v73
# Macaca mulatta	MMUL_1	v84
# Monodelphis domestica	BROADO5	v84
# Mus musculus	GRCm38	v84
# Oryctolagus cuniculus	OryCun2.0	v84
# Rattus norvegicus	Rnor_5.0	v79

# Homo_sapiens: wget https://ftp.ensembl.org/pub/release-73/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.73.dna.primary_assembly.fa.gz
# Macaca_mulatta: wget https://ftp.ensembl.org/pub/release-84/fasta/macaca_mulatta/dna/Macaca_mulatta.MMUL_1.dna.toplevel.fa.gz
# Gallus_gallus: wget https://ftp.ensembl.org/pub/release-84/fasta/gallus_gallus/dna/Gallus_gallus.Galgal4.dna.toplevel.fa.gz
# Monodelphis_domestica: https://ftp.ensembl.org/pub/release-84/fasta/monodelphis_domestica/dna/Monodelphis_domestica.BROADO5.dna.toplevel.fa.gz
# Mus_musculus: https://ftp.ensembl.org/pub/release-84/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
# Rattus_norvegicus: https://ftp.ensembl.org/pub/release-79/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_5.0.dna.toplevel.fa.gz
# Oryctolagus_cuniculus: wget https://ftp.ensembl.org/pub/release-84/fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa.gz


download_rna_seq_accession = {
    'human': 'E-MTAB-6814-unix-ftp.txt',
    'macaque': 'E-MTAB-6813-unix-ftp.txt',
    'mouse': 'E-MTAB-6798-unix-ftp.txt',
    'rat': 'E-MTAB-6811-unix-ftp.txt',
    'rabbit': 'E-MTAB-6782-unix-ftp.txt',
    'opossum': 'E-MTAB-6833-unix-ftp.txt',
    'chicken': 'E-MTAB-6769-unix-ftp.txt',
}

download_rna_seq = {
    'human': 'ftp://ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/814/E-MTAB-6814/Files',
    'macaque': 'ftp://ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/813/E-MTAB-6813/Files',
    'mouse': 'ftp://ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/798/E-MTAB-6798/Files',
    'rat': 'ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/811/E-MTAB-6811/Files',
    'rabbit': 'ftp://ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/782/E-MTAB-6782/Files',
    'opossum': 'ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/833/E-MTAB-6833/Files',
    'chicken': 'ftp://ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/769/E-MTAB-6769/Files',
}

download_gtf = {
    'human': 'https://ftp.ensembl.org/pub/release-73/gtf/homo_sapiens/Homo_sapiens.GRCh37.73.gtf.gz',
    'macaque': 'https://ftp.ensembl.org/pub/release-84/gtf/macaca_mulatta/Macaca_mulatta.MMUL_1.84.gtf.gz',
    'mouse': 'https://ftp.ensembl.org/pub/release-84/gtf/mus_musculus/Mus_musculus.GRCm38.84.gtf.gz',
    'rat': 'https://ftp.ensembl.org/pub/release-79/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_5.0.79.gtf.gz',
    'rabbit': 'https://ftp.ensembl.org/pub/release-84/gtf/oryctolagus_cuniculus/Oryctolagus_cuniculus.OryCun2.0.84.gtf.gz',
    'opossum': 'https://ftp.ensembl.org/pub/release-84/gtf/monodelphis_domestica/Monodelphis_domestica.BROADO5.84.gtf.gz',
    'chicken': 'https://ftp.ensembl.org/pub/release-84/gtf/gallus_gallus/Gallus_gallus.Galgal4.84.gtf.gz',
}

download_fasta = {
    'human': 'https://ftp.ensembl.org/pub/release-73/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.73.dna.primary_assembly.fa.gz',
    'macaque': 'https://ftp.ensembl.org/pub/release-84/fasta/macaca_mulatta/dna/Macaca_mulatta.MMUL_1.dna.toplevel.fa.gz',
    'mouse': 'https://ftp.ensembl.org/pub/release-84/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz',
    'rat': 'https://ftp.ensembl.org/pub/release-79/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_5.0.dna.toplevel.fa.gz',
    'rabbit': 'https://ftp.ensembl.org/pub/release-84/fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna.toplevel.fa.gz',
    'opossum': 'https://ftp.ensembl.org/pub/release-84/fasta/monodelphis_domestica/dna/Monodelphis_domestica.BROADO5.dna.toplevel.fa.gz',
    'chicken': 'https://ftp.ensembl.org/pub/release-84/fasta/gallus_gallus/dna/Gallus_gallus.Galgal4.dna.toplevel.fa.gz',
}

# get available timepoints for species and tissue
df_metadata = pd.read_csv(os.path.join(DATA_DIR, 'RNA_seq_files/metadata/df_annotation.tsv'), sep='\t')

grouped_df = df_metadata.groupby(['species', 'tissue', 'timepoint'])

df_metadata_stats = pd.DataFrame({
    'samples': grouped_df.apply(lambda df: sorted(set(df['library ID']))),
    'num_samples': grouped_df.apply(lambda df: len(set(df['library ID'])))
}).reset_index()

# Group by 'species' and 'tissue' and get the unique 'timepoint' values for each group
species_tissue_timepoints = df_metadata_stats.groupby(['species', 'tissue'])['timepoint'].unique().reset_index()

# Create a dictionary with ('species', 'tissue') as keys and 'timepoint' lists as values
species_tissue_dict = {
    (species, tissue): timepoints.tolist()
    for species, tissue, timepoints in species_tissue_timepoints.itertuples(index=False)
}


# create ID map for bam files of species, tissue, timepoint
from collections import defaultdict

# Define a nested defaultdict
# This will automatically create dictionaries at each level as needed
id_map = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

# Loop over the combinations of species, tissues, and timepoints
for s in species_kaesmann:
    for t in tissues_kaesmann:
        for tp in timepoints_kaesmann:

            # Select the relevant rows from the DataFrame
            _df = df_metadata[
                (df_metadata['species'] == s)
                & (df_metadata['tissue'] == t)
                & (df_metadata['timepoint'] == tp)
            ]

            # If there is data for this combination, store the library IDs
            if _df.shape[0] > 0:
                ids = sorted(set(_df['library ID']))
                ids = [str(x) for x in ids]
                
                # Store in the dictionary instead of a DataFrame
                id_map[s][t][tp] = ids