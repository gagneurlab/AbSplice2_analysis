##==================================OUTLIER GROUND TRUTH========================================= 
import pandas as pd

def tissue_DROP(wildcards):
    df_anno = pd.read_csv(config['DROP']['tissue_map'])
    tissue_map = dict(zip(df_anno['tissue'], df_anno['tissue_DROP']))
    return tissue_map[wildcards['tissue']]

def analysis_expression(wildcards):
    df_anno = pd.read_csv(config['DROP']['tissue_map'])
    tissue_map = dict(zip(df_anno['tissue'], df_anno['tissue_DROP']))
    return tissue_map[wildcards['tissue']] + '--' + config['DROP']['version']

def gene_expr_outlier(wildcards):
    tissue = tissue_DROP(wildcards)
    return config['DROP']['outlier']['expression']['working_dir'] + config['DROP']['version'] + f'/outrider/{tissue}/OUTRIDER_results.tsv'

# def fraser1_path_orig(wildcards):
#     df_anno = pd.read_csv(config['DROP']['tissue_map'])
#     tissue_map = dict(zip(df_anno['tissue'], df_anno['tissue_DROP']))
#     return config['DROP']['outlier_ground_truth_all'] + tissue_map[wildcards['tissue']] + '/results_per_junction_fdr.tsv'

def fraser_path(wildcards):
    df_anno = pd.read_csv(config['DROP']['tissue_map'])
    tissue_map = dict(zip(df_anno['tissue'], df_anno['tissue_DROP']))
    if wildcards.fraser_version == 'fraser1':
        path_to_add = config['DROP']['fraser1_results']
    elif wildcards.fraser_version == 'fraser2':
        path_to_add = config['DROP']['fraser2_results']
    path = os.path.join(
        path_to_add,
        config['DROP']['version'],
        'fraser',
        tissue_map[wildcards['tissue']]
    )
    return path

# def fraser2_path(wildcards):
#     df_anno = pd.read_csv(config['DROP']['tissue_map'])
#     tissue_map = dict(zip(df_anno['tissue'], df_anno['tissue_DROP']))
#     path = os.path.join(
#         config['DROP']['fraser2_results'],
#         config['DROP']['version'],
#         'fraser',
#         tissue_map[wildcards['tissue']]
#     )
#     return path
        

# # Unified Outlier Processing Rule
# rule outlier_processing:
#     input:
#         results=lambda wildcards: fraser1_path(wildcards) if wildcards.fraser_version == 'fraser1' else fraser2_path(wildcards)
#     output:
#         junction_level = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['fraser']['junction_level'],
#         gene_level = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['fraser']['gene_level'],
#     conda:
#         lambda wildcards: '../../../../../envs/environment_FRASER.yaml' if wildcards.fraser_version == 'fraser1' else '../../../../../envs/environment_FRASER2.yaml'
#     resources:
#         mem_mb=lambda wildcards, attempt: attempt * 265000,
#         threads=4
#     script:
#         lambda wildcards: "../fraser_outlier_ground_truth_unfiltered.R" if wildcards.fraser_version == 'fraser1' else "../fraser2_gene_junction_results.py"


# rule fraser1_orig:
#     input:
#         results = fraser1_path_orig
#     params:
#         tissue = '{tissue}',
#         # delta_psi_cutoff = config['filtering_params']['outliers']['delta_psi_cutoff'],
#         # outlier_type = config['filtering_params']['outliers']['outlier_type'],
#     conda:
#         '../../../../../envs/environment_FRASER.yaml'
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 265000,
#         threads = 4
#     output:
#         junction_level = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['fraser']['junction_level'],
#         gene_level = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['fraser']['gene_level'],
#     script:
#         "../fraser_outlier_ground_truth_unfiltered.R"


# if 'fraser2' in config['fraser_version']:
rule fraser:
    input:
        results = fraser_path # TODO: depending on DROP version the gene level file is called results_gene_all.tsv or results.tsv. results.tsv is already filtered for significant
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 265000,
        threads = 4
    # params:
    #     delta_psi_cutoff = config['filtering_params']['outliers']['delta_psi_cutoff'],
    output:
        junction_level = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['fraser']['junction_level'],
        gene_level = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['fraser']['gene_level'],
    script:
        "../fraser2_gene_junction_results.py"


rule outliers_gene_level:
    input:
        gene_level = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['fraser']['gene_level'],
        sample_map = config['DROP']['sample_annotation'],
        gene_map = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
        coding_genes = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['coding_genes'],
    params:
        tissue_DROP = tissue_DROP,
        # padjustGene_cutoff = config['filtering_params']['outliers']['padjustGene_cutoff'],
        key_assay = config['filtering_params']['count_table']['key_assay'],
        value_assay = config['filtering_params']['count_table']['value_assay'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['qual_filtered']['gene_level'],
    script:
        "../outliers_gene_level.py"
        
        
rule outliers_junction_level_signif:
    input:
        junction_level = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['fraser']['junction_level'],
        sample_map = config['DROP']['sample_annotation'],
        count_table_updated = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['count_table']['updated'],
        gene_map = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
        coding_genes = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['coding_genes'],
    params:
        tissue_DROP = tissue_DROP,
        # padjust_junction_cutoff = config['filtering_params']['outliers']['padjust_junction_cutoff'],
        # totalCounts_cutoff = config['filtering_params']['outliers']['totalCounts_cutoff'],
        # delta_psi_cutoff = config['filtering_params']['outliers']['delta_psi_cutoff'],
        num_junction_outliers_per_samples = config['filtering_params']['outliers']['num_junction_outliers_per_samples'],
        key_assay = config['filtering_params']['count_table']['key_assay'],
        value_assay = config['filtering_params']['count_table']['value_assay'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        junction_signif = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['qual_filtered']['junction_level_signif'],
    script:
        "../outliers_junction_level_signif.py"

        
rule outliers_combine_gene_junction_level_signif:
    input:
        junction = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['qual_filtered']['junction_level_signif'],
        gene = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['qual_filtered']['gene_level'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 1
    output:
        gene_junction = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['gene_junction_signif'],
    script:
        "../combine_gene_junction.py"
                 
            
# outliers with distance to variants
rule outlier_var_nearest_variant_centric:
    input:
        count_table_annotation = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['count_table']['with_annotation'],
        outliers_signif = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['gene_junction_signif'],
        rare_vars = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    output:
        var_junc_nearest = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['parts']
    script:
        "../outlier_var_nearest_variant_centric.py.py"

        
rule outlier_var_nearest_variant_centric_rare_var_dist:
    input:
        outliers = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['parts'],
    output:
        outliers_subset = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['parts_rare_var_dist'],
    run:
        df = pd.read_csv(input.outliers)
        df = df[
            df['abs_Distance'] <= 250]
        df.to_csv(output.outliers_subset, index=False)
        
        
rule outlier_var_nearest_variant_centric_rare_var_dist_gene_level:
    input:
        outliers_variant_level = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['parts_rare_var_dist'],
    output:
        outliers_gene_level = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['parts_rare_var_dist_gene_level'],
    run:
        df = pd.read_csv(input.outliers_variant_level)
        df = df[['gene_id', 'sample']].drop_duplicates()
        df.to_csv(output.outliers_gene_level, index=False)
        
        
rule outlier_var_nearest_variant_centric_rare_var_dist_variant_level:
    input:
        outliers_variant_level = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['parts_rare_var_dist'],
    output:
        outliers_variant_level = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['parts_rare_var_dist_variant_level'],
    run:
        df = pd.read_csv(input.outliers_variant_level)
        df = df[['variant', 'gene_id', 'sample']].drop_duplicates()
        df.to_csv(output.outliers_variant_level, index=False)
        
        
rule outliers_tissue_cat_pval:
    input:
        outliers_pval = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['parts_rare_var_dist'],
        outliers_variant_level = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['parts_rare_var_dist_variant_level'],
    output:
        outliers_pval = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['tissue_cat_pval'],
    script:
        "../outliers_tissue_cat_pval.py"


rule outliers_tissue_cat_pval_gene_level:
    input:
        outliers_pval_variant_level = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['tissue_cat_pval'],
    output:
        outliers_pval_gene_level = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['tissue_cat_pval_gene_level'],
    script:
        "../outliers_tissue_cat_pval_gene_level.py"
        
import re
valid_combinations = config['filtering_params']['outliers']['combinations']

def filter_wildcards_outliers(path_list, valid_combinations):
    """Filters paths to include only valid combinations of delta_psi_cutoff, padjust_cutoff, and totalCounts_cutoff."""
    correct_paths = []
    
    for path in path_list:
        match = re.search(
            r"deltaPsi=(?P<delta_psi_cutoff>[\d.]+)/pval=(?P<padjust_cutoff>[\d.]+)/totalCounts=(?P<totalCounts_cutoff>[\d.]+)", 
            path
        )
        
        if match:
            delta_psi_cutoff = float(match.group("delta_psi_cutoff"))
            padjust_cutoff = float(match.group("padjust_cutoff"))
            totalCounts_cutoff = float(match.group("totalCounts_cutoff"))

            # Check if the extracted wildcards match any valid combination
            if any(
                c["delta_psi_cutoff"] == delta_psi_cutoff and
                c["padjust_cutoff"] == padjust_cutoff and
                c["totalCounts_cutoff"] == totalCounts_cutoff
                for c in valid_combinations
            ):
                correct_paths.append(path)
    
    return correct_paths


rule all_outlier_ground_truth:
    input:
        filter_wildcards_outliers(expand(
            rules.fraser.output, 
            tissue=config['tissue_target'], 
            fraser_version=config['fraser_version'], 
            delta_psi_cutoff=[c['delta_psi_cutoff'] for c in valid_combinations],
            padjust_cutoff=[c['padjust_cutoff'] for c in valid_combinations],
            totalCounts_cutoff=[c['totalCounts_cutoff'] for c in valid_combinations],
        ), valid_combinations),
        filter_wildcards_outliers(expand(
            rules.outlier_var_nearest_variant_centric_rare_var_dist_variant_level.output, 
            vcf_id=wildcard_vcf_id,
            tissue=config['tissue_target'], 
            fraser_version=config['fraser_version'], 
            delta_psi_cutoff=[c['delta_psi_cutoff'] for c in valid_combinations],
            padjust_cutoff=[c['padjust_cutoff'] for c in valid_combinations],
            totalCounts_cutoff=[c['totalCounts_cutoff'] for c in valid_combinations],
        ), valid_combinations)