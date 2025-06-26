import pandas as pd
from splicemap.splice_map import SpliceMap

def valid_splicemaps(input_files, output_file):
    splicemaps = []
    for filename in input_files:
        try:
            SpliceMap.read_csv(filename)
            splicemaps.append(filename)
        except:
            pass
    splicemaps = [str(x) for x in splicemaps]

    # Open a file in write mode
    with open(output_file, 'w') as file:
        for item in splicemaps:
            file.write(f"{item}\n")
            
            
valid_splicemaps(snakemake.input['splicemaps5'], snakemake.output['valid_splicemaps5'])
valid_splicemaps(snakemake.input['splicemaps3'], snakemake.output['valid_splicemaps3'])

