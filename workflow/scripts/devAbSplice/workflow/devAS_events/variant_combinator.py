from kipoiseq.extractors import VariantCombinator

VariantCombinator(
    snakemake.input['fasta'], snakemake.input['bed']
).to_vcf(snakemake.output['vcf'])
