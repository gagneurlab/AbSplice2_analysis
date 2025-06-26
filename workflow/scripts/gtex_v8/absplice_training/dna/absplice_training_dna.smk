import os
import yaml
import pandas as pd
from tqdm import tqdm
from absplice_scripts.data.gtex_v8 import tissue_map, tissue_map_gene_expr

ABSPLICE_INPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['combined_input']

#------------------------------------------- DNA --------------------------------------
rule absplice_dna_training_input:
    input:
        universe = OUTPUT_DIR_BENCHMARK + config_static['benchmark']['universe'],
        outliers = OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['parts_rare_var_dist_variant_level'],
        absplice_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
    params:
        tissue = '{tissue}',
        vcf_id = '{vcf_id}',
    output:
        absplice_input = OUTPUT_DIR_SPLICING + config['absplice_training']['training_input']['dna']
    script:
        "../absplice_combine_training_input.py"
        
        
rule absplice_dna_training_input_subset:
    input:
        absplice_input = OUTPUT_DIR_SPLICING + config['absplice_training']['training_input']['dna']
    params:
        tissue = '{tissue}',
        vcf_id = '{vcf_id}',
    output:
        absplice_input = OUTPUT_DIR_SPLICING + config['absplice_training']['training_input_subset']['dna']
    script:
        "./training_input_subset_dna.py"
        
        
rule absplice_dna_input_complete:
    input:
        absplice_input = expand(OUTPUT_DIR_SPLICING + config['absplice_training']['training_input_subset']['dna'],
                                tissue=config['tissues'], vcf_id=wildcard_vcf_id,
                                fraser_version='{fraser_version}',
                                delta_psi_cutoff= '{delta_psi_cutoff}',
                                padjust_cutoff='{padjust_cutoff}',
                                totalCounts_cutoff='{totalCounts_cutoff}'
                                )
    output:
        absplice_input = OUTPUT_DIR_SPLICING + config['absplice_training']['training_input_complete']['dna']
    run:
        df = pd.concat([pd.read_csv(i) for i in tqdm(input.absplice_input)])
        df.to_csv(output.absplice_input, index=False)
    
    
rule absplice_dna_training_5fold:
    input:
        absplice_input = OUTPUT_DIR_SPLICING + config['absplice_training']['training_input_complete']['dna'],
    params:
        classifier = '{classifier}',
        feature_string = '{feature_string}',
        abs_features = '{abs_features}',
        median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
        gene_tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
    output:
        absplice_5fold_crossval = OUTPUT_DIR_SPLICING + config['absplice_training']['preds']['dna']['5_fold_crossval'],
    script:
        "../absplice_training_5fold.py"
        
        
rule absplice_dna_training_5fold_gene_level:
    input:
        pred_absplice_dna = OUTPUT_DIR_SPLICING + config['absplice_training']['preds']['dna']['5_fold_crossval'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_absplice_dna_gene = OUTPUT_DIR_SPLICING + config['absplice_training']['preds']['dna']['gene_level'],
    script:
        "../gene_level_absplice_dna.py"
        
        
rule absplice_dna_training_5fold_variant_level:
    input:
        pred_absplice_dna = OUTPUT_DIR_SPLICING + config['absplice_training']['preds']['dna']['5_fold_crossval'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1
    output:
        result_absplice_dna_variant = OUTPUT_DIR_SPLICING + config['absplice_training']['preds']['dna']['variant_level'],
    script:
        "../variant_level_absplice_dna.py"
        
        
rule absplice_dna_training_whole_GTEx:
    input:
        absplice_input = OUTPUT_DIR_SPLICING + config['absplice_training']['training_input_complete']['dna'],
    params:
        classifier = '{classifier}',
        feature_string = '{feature_string}',
        abs_features = '{abs_features}',
        median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
        gene_tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
    output:
        absplice_whole_GTEx = OUTPUT_DIR_SPLICING + config['absplice_training']['preds']['dna']['whole_GTEx'],
    script:
        "../absplice_training_whole_GTEx.py"
        

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

        
rule all_absplice_training_dna:
    input:
        filter_wildcards_outliers(expand(rules.absplice_dna_training_5fold_gene_level.output,
               classifier=config['absplice_params']['classifier'],
               abs_features=config['absplice_params']['abs_features'],
               feature_string=config['absplice_params']['feature_string_dna'],
               fraser_version=config['fraser_version'], 
                delta_psi_cutoff=[c['delta_psi_cutoff'] for c in valid_combinations],
                padjust_cutoff=[c['padjust_cutoff'] for c in valid_combinations],
                totalCounts_cutoff=[c['totalCounts_cutoff'] for c in valid_combinations],
            ), valid_combinations),
        filter_wildcards_outliers(expand(rules.absplice_dna_training_5fold_variant_level.output,
               classifier=config['absplice_params']['classifier'],
               abs_features=config['absplice_params']['abs_features'],
               feature_string=config['absplice_params']['feature_string_dna'],
               fraser_version=config['fraser_version'], 
                delta_psi_cutoff=[c['delta_psi_cutoff'] for c in valid_combinations],
                padjust_cutoff=[c['padjust_cutoff'] for c in valid_combinations],
                totalCounts_cutoff=[c['totalCounts_cutoff'] for c in valid_combinations],
            ), valid_combinations),
        filter_wildcards_outliers(expand(rules.absplice_dna_training_whole_GTEx.output,
               classifier=config['absplice_params']['classifier'],
               abs_features=config['absplice_params']['abs_features'],
               feature_string=config['absplice_params']['feature_string_dna'],
               fraser_version=config['fraser_version'], 
                delta_psi_cutoff=[c['delta_psi_cutoff'] for c in valid_combinations],
                padjust_cutoff=[c['padjust_cutoff'] for c in valid_combinations],
                totalCounts_cutoff=[c['totalCounts_cutoff'] for c in valid_combinations],
            ), valid_combinations),


del valid_combinations