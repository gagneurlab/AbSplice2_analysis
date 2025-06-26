from absplice import SplicingOutlierResult

splicing_result = SplicingOutlierResult(
        df_mmsplice=snakemake.input['mmsplice_splicemap'], 
        df_spliceai=snakemake.input['spliceai'],
        gene_map=snakemake.params['gene_map'],
    )
print('Loaded input')
splicing_result.predict_absplice_dna(
    pickle_file=snakemake.params['absplice_dna_model'], 
    extra_info=False)
print('Saving results')
splicing_result._absplice_dna.reset_index().to_parquet(
    snakemake.output['absplice_dna'], 
    partition_cols=['tissue'],
    index=False)