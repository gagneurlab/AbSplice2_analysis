import pandas as pd
import pyranges as pr
from splicemap import SpliceCountTable as CountTable
from splicemap.splice_map import SpliceMap
from tqdm import tqdm
tqdm.pandas()

s = 'human'
t = snakemake.wildcards['tissue']
tp = snakemake.wildcards['timepoint']

fasta_file = snakemake.input['fasta_file']
gtf_file = snakemake.input['gtf_file']

# read in count table
ct = CountTable.read_csv(snakemake.input['ct'], name=f'{s}_{t}_{tp}')

if ct.df.shape[0] > 0:
    # # remove chr annotation
    # if snakemake.params['remove_chr']:
    #     ct.df['Chromosome'] = ct.df['Chromosome'].str.replace('chr', '')
    #     ct.df = ct.df.reset_index()
    #     ct.df['junctions'] = ct.df['junctions'].str.replace('chr', '')
    #     ct.df = ct.df.set_index('junctions')
        
    # infer strand
    ct_no_strand = ct.df.copy()
    ct_no_strand['Strand'] = '.'
    ct = CountTable(ct_no_strand, name=f'{s}_{t}_{tp}')
    ct.infer_strand(fasta_file, progress=True)

    # drop duplicates
    def process_group(df):
        first_part = df[['Chromosome', 'Start', 'End', 'Strand']].iloc[0]
        second_part = df.iloc[:, 4:].sum()
        return pd.concat([first_part, second_part])

    _df = ct.df
    _df = _df.groupby('junctions').progress_apply(process_group)
    # _df = _df.groupby('junctions').progress_apply(
    #     lambda df: df[['Chromosome', 'Start', 'End', 'Strand']].iloc[0].append(df.iloc[:, 4:].sum()))
    ct = CountTable(_df, name=f'{s}_{t}_{tp}')

    # median_n_cutoff
    ct_psi5 = ct.event5_median_filter(cutoff=snakemake.params['median_n_cutoff'])
    ct_psi3 = ct.event3_median_filter(cutoff=snakemake.params['median_n_cutoff'])

    # quantile filter (in splicemap np.percentile behaves strange with read counts, e.g. if 0,0,0,2 passes, if 0,0,0,1 does not pass)
    def _quantile_filter(counts, quantile=95, min_read=1):
        # Calculate the percentage of samples with reads >= min_read for each junction
        percentage_filter = (counts >= min_read).sum(axis=1) / counts.shape[1] * 100
        # Filter junctions that are observed in at least the given quantile of samples
        filtered_junctions = percentage_filter >= (100 - quantile)
        return filtered_junctions
    # do not need to filter, because already filtered in subset_count_table.py (at least 1 sample with >= 1 read)
    ct_psi5.df = ct_psi5.df[_quantile_filter(ct_psi5.df.iloc[:, 4:])]
    ct_psi3.df = ct_psi3.df[_quantile_filter(ct_psi3.df.iloc[:, 4:])]

    # infer annotation
    ct_psi5.infer_annotation(gtf_file)
    ct_psi3.infer_annotation(gtf_file)

    # ref psi
    df_psi5 = ct_psi5.ref_psi5(method='k/n')
    df_psi3 = ct_psi3.ref_psi3(method='k/n')

    # save
    df_psi5.to_csv(snakemake.output['splicemap_psi5'])
    df_psi3.to_csv(snakemake.output['splicemap_psi3'])

else:
    # save empty dataframe
    cols = ['junctions', 'gene_id', 'Chromosome', 'Start', 'End', 'Strand',
       'splice_site', 'events', 'ref_psi', 'k', 'n', 'median_n', 'median_k',
       'mean_n', 'mean_k', 'gene_name', 'gene_type', 'novel_junction',
       'weak_site_donor', 'weak_site_acceptor', 'transcript_id', 'gene_tpm']
    df_empty = pd.DataFrame(columns=cols)
    SpliceMap(
        df_empty,
        f'{s}_{t}_{tp}'
    ).to_csv(snakemake.output['splicemap_psi5'])
    SpliceMap(
        df_empty,
        f'{s}_{t}_{tp}'
    ).to_csv(snakemake.output['splicemap_psi3'])