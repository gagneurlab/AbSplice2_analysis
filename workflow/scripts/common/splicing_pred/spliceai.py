# import pyarrow
from spliceai_rocksdb.spliceAI import SpliceAI

# genome_map = {
#     'hg19': 'grch37',
#     'hg38': 'grch38'
# }

# if snakemake.params['lookup_only']:
#     model = SpliceAI(db_path=snakemake.input['db'])
# else:
#     model = SpliceAI(fasta=snakemake.input['fasta'],
#                      annotation=genome_map[snakemake.params['genome']],
#                      db_path=snakemake.input['db'])

if snakemake.params['lookup_only']:
    model = SpliceAI(db_path=snakemake.params['db_path'])
else:
    model = SpliceAI(snakemake.input['fasta'],
                     annotation=snakemake.params['genome'],
                     db_path=snakemake.params['db_path'])

model.predict_save(snakemake.input['vcf'],
                   snakemake.output['result'])