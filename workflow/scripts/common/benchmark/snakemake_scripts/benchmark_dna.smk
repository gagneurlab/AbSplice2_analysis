from absplice_scripts.utils.model_utils import *

OUTPUT_DIR_TISSUE_SUBSET = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['tissue_subset'] + config_static['splicing_pred']['levels']['gene_level']
OUTPUT_DIR_GENE = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['gene_level']

# -----------------------------------GTEx SpliceMaps-------------------------------------
# def get_model_input(model_name, wildcards):
#     """Returns the path for the specified model if it is present in config['models'], otherwise None."""
#     output_paths = {
#         "spliceai": OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['spliceai'],
#         "spliceai_splicemap": OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['spliceai_splicemap']['all'],
#         "spliceai_splicemap_ref_psi": OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['all'],
#         "mmsplice_splicemap": OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['mmsplice_splicemap'],
#         "mmsplice_splicemap_ref_psi": OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['mmsplice_splicemap_ref_psi'],
#         "absplice_dna": OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['absplice']['dna'],
#     }
#     return output_paths.get(model_name) if model_name in config['models'] else None

# rule benchmark_combine_dna:
#     input:
#         universe=OUTPUT_DIR_BENCHMARK + config_static['benchmark']['universe'],
#         outliers=OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['parts_rare_var_dist_gene_level'],
#         spliceai=lambda wildcards: get_model_input('spliceai', wildcards),
#         spliceai_splicemap=lambda wildcards: get_model_input('spliceai_splicemap', wildcards),
#         spliceai_splicemap_ref_psi=lambda wildcards: get_model_input('spliceai_splicemap_ref_psi', wildcards),
#         mmsplice_splicemap=lambda wildcards: get_model_input('mmsplice_splicemap', wildcards),
#         mmsplice_splicemap_ref_psi=lambda wildcards: get_model_input('mmsplice_splicemap_ref_psi', wildcards),
#         absplice_dna=lambda wildcards: get_model_input('absplice_dna', wildcards),
#     resources:
#         mem_mb=lambda wildcards, attempt: attempt * 64000,
#     params:
#         pred_tools=pred_tools_dna,
#         tissue_specific_tools=tissue_specific_tools,
#         tissue_pred='{tissue_pred}',
#         tissue='{tissue}',
#         cols_spliceai=['variant', 'delta_score'],
#         cols_mmsplice=['variant', 'delta_logit_psi'],
#         cols_mmsplice_splicemap=['variant', 'delta_logit_psi', 'delta_psi', 'median_n', 'ref_psi', 'tissue'],
#         cols_mmsplice_splicemap_ref_psi=['variant', 'delta_logit_psi', 'delta_psi', 'median_n', 'ref_psi', 'tissue'],
#         cols_absplice_dna=['variant', 'AbSplice_DNA', 'tissue', 'delta_logit_psi', 'delta_psi', 'delta_score', 'splice_site_is_expressed'],
#     output:
#         combined_benchmark=OUTPUT_DIR_BENCHMARK + config_static['benchmark']['combined_benchmark']['dna'],
#     script:
#         "../benchmark_combine_dna.py"

# cols_absplice_custom_new = []
cols_absplice_custom_new = ['tissue']
for model_name, model_info in config_static['absplice_models_new'].items():
    _cols = [
        f'variant_{model_name}', model_name, 
        *[f'{x}_{model_name}' for x in model_info['model_features']]
    ]
    cols_absplice_custom_new.extend(_cols)
    model_dict_dna[model_name] = f'{model_name}_absplice_custom_new'


rule benchmark_combine_dna:
    input:
        universe = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['universe'],
        outliers = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['parts_rare_var_dist_gene_level'],
        
        spliceai = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['spliceai'],
        # mmsplice_splicemap = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['mmsplice_splicemap'],
        absplice_dna = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['absplice']['dna'],
        # spliceai_splicemap = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['spliceai_splicemap']['all'],
        # spliceai_splicemap_ref_psi = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['all'],
        # mmsplice_splicemap_ref_psi = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['mmsplice_splicemap_ref_psi'],
        pangolin = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['pangolin'],
        absplice_custom_new = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['absplice_custom_new']['dna'],

        # mmsplice = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['mmsplice'],
        # mtsplice = OUTPUT_DIR_TISSUE_SUBSET + config_static['splicing_pred']['models']['mtsplice'],
        # cadd_splice = OUTPUT_DIR_GENE + config_static['splicing_pred']['models']['cadd_splice'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    params:
        pred_tools = pred_tools_dna,
        tissue_specific_tools = tissue_specific_tools,
        tissue_pred = '{tissue_pred}',
        tissue = '{tissue}',
        cols_absplice_custom_new = cols_absplice_custom_new,
        # cols_absplice_custom_new = [
        #     # 'variant_model_A', 'model_A', 
        #     'variant_model_B', 'model_B',
        #      'tissue'],
        cols_spliceai = [
            'variant', 'delta_score'],
        cols_pangolin = [
            'variant', 'Pangolin_max_score'],
        cols_mmsplice = [
            'variant', 'delta_logit_psi'],
        cols_mmsplice_splicemap = [
            'variant', 'delta_logit_psi', 'delta_psi', 'median_n', 'ref_psi', 'tissue'],
        cols_mmsplice_splicemap_ref_psi = [
            'variant', 'delta_logit_psi', 'delta_psi', 'median_n', 'ref_psi', 'tissue'],
        cols_absplice_dna = [
            'variant', 'AbSplice_DNA', 'tissue', 
            'delta_logit_psi', 'delta_psi', 'delta_score', 'splice_site_is_expressed'],
        cols_mtsplice = [
            'variant', 'delta_logit_psi_mtsplice', 'tissue'
        ],
        cols_cadd_splice = [
            'variant', 'PHRED'
        ],
        cols_spliceai_splicemap = [
            'junctions', 'delta_score', 'psi5_median_n', 'psi3_median_n', 'psi5_ref_psi', 'psi3_ref_psi', 'tissue'
        ],
        cols_spliceai_splicemap_ref_psi = [
            'junctions', 'delta_score', 'psi5_median_n', 'psi3_median_n', 'psi5_ref_psi', 'psi3_ref_psi', 'tissue'
        ]
    output:
        combined_benchmark = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['combined_benchmark']['dna'],
    script:
        "../benchmark_combine_dna.py"
        
        
rule performance_dna_tissue_subset:
    input:
        benchmark = expand(OUTPUT_DIR_BENCHMARK + config_static['benchmark']['combined_benchmark']['dna'],
                           vcf_id=wildcard_vcf_id, tissue='{tissue}', tissue_pred='{tissue_pred}', fraser_version='{fraser_version}',
                           delta_psi_cutoff='{delta_psi_cutoff}', padjust_cutoff='{padjust_cutoff}', totalCounts_cutoff='{totalCounts_cutoff}'),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
    params:
        outlier_column = 'outlier',
        model_dict = model_dict_dna,
        median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
        gene_tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
    output:
        df_performance = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['performance']['dna']['df'],
        aps_performance = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['performance']['dna']['aps']
    script:
        "../performance_dna.py"
        
        
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


rule all_benchmark_dna:
    input:
        filter_wildcards_outliers(expand(rules.performance_dna_tissue_subset.output,
               tissue=config['tissue_target'],
            #    tissue_pred=config['tissues'],
               tissue_pred=config['tissues_subset'],
               fraser_version=config['fraser_version'],
                #    delta_psi_cutoff=config['filtering_params']['outliers']['delta_psi_cutoff'],
                #     padjust_cutoff=config['filtering_params']['outliers']['padjust_junction_cutoff'],
                #     totalCounts_cutoff=config['filtering_params']['outliers']['totalCounts_cutoff'],
                delta_psi_cutoff=[c['delta_psi_cutoff'] for c in valid_combinations],
                padjust_cutoff=[c['padjust_cutoff'] for c in valid_combinations],
                totalCounts_cutoff=[c['totalCounts_cutoff'] for c in valid_combinations],
              ), valid_combinations)
        # expand(rules.performance_dna_tissue_subset.output,
        #        tissue=config['tissue_target'],
        #        tissue_pred=config['tissues_subset'] 
        #       ),
        
        
