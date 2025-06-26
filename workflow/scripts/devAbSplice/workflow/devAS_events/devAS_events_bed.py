import pandas as pd
import pyranges as pr

def convert_strand(strand):
    if strand == -1:
        return '-'
    elif strand == 1:
        return '+'
    else:
        print('missin')
        raise KeyError()
    
def create_bed(_df, window_size = 100, chrom='chr1'):
    df = _df[['chr_id', 'start', 'stop', 'strand']]
    df = df.rename(columns={
        'chr_id': 'Chromosome',
        'start': 'Exon_start',
        'stop': 'Exon_end',
        'strand': 'Strand'
    })
    df['Chromosome'] = df['Chromosome'].apply(lambda df: 'chr' + str(df))
    df['Strand'] = df.apply(lambda df: convert_strand(df['Strand']), axis=1)
    
    # concat start and end of exon
    df = pd.concat([
        df[['Chromosome', 'Exon_start', 'Strand']].rename(columns={'Exon_start': 'pos'}),
        df[['Chromosome', 'Exon_end', 'Strand']].rename(columns={'Exon_end': 'pos'})
    ])
    
    df['Start'] = df['pos'] - window_size
    df['End'] = df['pos'] + window_size
    
    df = df[['Chromosome', 'Start', 'End', 'Strand', 'pos']]
    df = df[
        df['Chromosome'] == chrom
    ]
    assert len(set(df['Chromosome'])) <= 1
    pr_region = pr.PyRanges(df)
    return pr_region

# create bed
df = pd.read_csv(snakemake.input['devAS_events'])
pr_region = create_bed(df, window_size = 100, chrom=snakemake.wildcards['chrom'])
pr_region.to_bed(snakemake.output['devAS_events_bed'])

