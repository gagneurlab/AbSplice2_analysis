# DIR_GTEX = '../../../../AbSplice_analysis/workflow/data/results/gtex_v8/junction_annotation/'
DIR_GTEX = '../../data/results/gtex_v8/junction_annotation/'
DIR_GVEX = '../../data/results/gvex/junction_annotation/'


rule copy_gene_tpm_from_gtex:
    input:
        gene_tpm = os.path.join(DIR_GTEX, 'gene_tpm.csv'),
    output:
        gene_tpm = os.path.join(DIR_GVEX, 'gene_tpm.csv'),
    shell:
        """
        cp {input.gene_tpm} {output.gene_tpm}
        """