OUTPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['raw']
ABSPLICE_INPUT_DIR = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['combined_input']

##==================================SPLICING PREDICTIONS (BASELINE MODELS)=========================================    
# MMSplice
rule splicing_pred_mmsplice:
    input:
        vcf = config['vcf'],
        fasta = config['fasta'],
        gtf = config['gtf'],
    params:
        genome = config['genome'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1
    output:
        result = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice'])
        # result = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice']
    script:
        "../mmsplice_exon.py"
        
# MTSplice
rule splicing_pred_mtsplice:
    input:
        vcf = config['vcf'],
        fasta = config['fasta'],
        gtf = config['gtf'],
    params:
        genome = config['genome'],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1
    output:
        result = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['mtsplice'])
        # result = OUTPUT_DIR + config_static['splicing_pred']['models']['mtsplice']
    script:
        "../mtsplice.py"

# # SpliceAI
# rule splicing_pred_spliceai:
#     input:
#         vcf = config['vcf'],
#         db = config_precomputed['spliceai']['db'].format(
#             genome=config['genome']),
#         fasta = config['fasta'],
#     conda:
#         '../../../../../envs/environment_spliceai_rocksdb.yaml'
#     params:
#         lookup_only = config['spliceai']['lookup_only'],
#         genome = config['genome'],
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 16000,
#         threads = 1,
#         gpu = 1,
#     output:
#         result = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'])
#         # result = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai']
#     script:
#         "../spliceai.py"
#     # conda: 'spliceai-rocksdb_test'


genome_mapper = {
    'hg38': 'grch38',
    'hg19': 'grch37',
}

def dict_path(wildcards):
    paths = {}
    genome =config['genome']
    for chr in config['chroms']:
        paths[chr] = config_precomputed['spliceai']['db'].format(genome=genome_mapper[genome], chrom=chr)
    return paths

rule splicing_pred_spliceai:
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1,
        # gpu = 1,
    input:
        vcf = config['vcf'],
        fasta = config['fasta'],
        db = expand(config_precomputed['spliceai']['db'],
                genome=genome_mapper[config['genome']], chrom=config['chroms']),
    params:
        db_path = dict_path,
        lookup_only = config['spliceai']['lookup_only'],
        genome = genome_mapper[config['genome']],
    conda:
        '../../../../../envs/environment_spliceai_rocksdb.yaml'
    output:
        result = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'])
    script:
        "../spliceai.py"


# Pangolin
rule splicing_pred_pangolin:
    input:
        vcf = config['vcf'],
        fasta = config['fasta'],
        anno = config['anno_pangolin']
    # conda: 'link_to_yaml' # TODO: need to implement
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1,
        gpu = 1,
    output:
        result = OUTPUT_DIR + config_static['splicing_pred']['models']['pangolin_raw']
    shell:
        'pangolin -m True -d 50 {input.vcf} {input.fasta} {input.anno} {output.result}'
     
        
# CADD-Splice        
genome_mapper = {
    'hg38': 'GRCh38',
    'hg19': 'GRCh37',
}

rule splicing_pred_cadd_splice:
    input:
        vcf = config['vcf'],
    conda: '/opt/modules/i12g/anaconda/envs/cadd_splice'
    params:
        genome = genome_mapper[config['genome']],
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        threads = 1
    output:
        result = OUTPUT_DIR + config_static['splicing_pred']['models']['cadd_splice']
    shell:
        '''/data/nasif12/home_if12/wagnern/Projects/CADD-scripts/CADD.sh -a -g "{params.genome}" -o {output.result} {input.vcf}'''
    

##==================================SPLICING PREDICTIONS (WITH SPLICEMAP)=========================================    
# MMSplice + SpliceMap
rule splicing_pred_mmsplice_splicemap:
    input:
        vcf = config['vcf'],
        fasta = config['fasta'],
        splicemap_5 = expand(config_static['junction_annotation']['splicemap']['psi5'],
                             genome=config['genome'], tissue=config['tissues']),
        splicemap_3 = expand(config_static['junction_annotation']['splicemap']['psi3'],
                             genome=config['genome'], tissue=config['tissues']),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 64000,
        threads = 4
    output:
        result = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap'])
        # result = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap']
    script:
        "../mmsplice_splicemap.py"
        
# SpliceAI + SpliceMap (single tissue)
rule splicing_pred_spliceai_splicemap_tissue:
    input:
        spliceai = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'],
        splicemap_5 = expand(config_static['junction_annotation']['splicemap']['psi5'],
                             genome=config['genome'], tissue='{tissue}'),
        splicemap_3 = expand(config_static['junction_annotation']['splicemap']['psi3'],
                             genome=config['genome'], tissue='{tissue}'),
        gene_map = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
    params:
        tissue = '{tissue}'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1,
    output:
        spliceai_splicemap = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['tissue'])
        # spliceai_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['tissue']
    script:
        "../spliceai_splicemap.py"

# SpliceAI + SpliceMap (all tissues)      
rule splicing_pred_spliceai_splicemap_all:
    input:
        df_tissue = expand(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['tissue'],
                           vcf_id='{vcf_id}', tissue=config['tissues']),
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1,
    output:
        df_all = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['all'])
        # df_all = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['all']
    run:
        import pandas as pd
        from tqdm import tqdm
        df = pd.concat([pd.read_parquet(i) for i in tqdm(input.df_tissue)])
        df.to_parquet(output.df_all, index=False, partition_cols=['tissue'])
        
# SpliceAI + SpliceMap + PSI_ref (single tissue)    
rule splicing_pred_spliceai_splicemap_ref_psi_tissue:
    input:
        spliceai_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap']['tissue']
    params:
        tissue = '{tissue}',
        ref_psi_tolerance = 0.05
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        threads = 1,
    output:
        spliceai_splicemap_ref_psi = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['tissue'])
        # spliceai_splicemap_ref_psi = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['tissue']
    script:
        "../spliceai_splicemap_ref_psi.py"
        
# SpliceAI + SpliceMap + PSI_ref (all tissues)         
rule splicing_pred_spliceai_splicemap_ref_psi_all:
    input:
        df_tissue = expand(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['tissue'],
                           vcf_id='{vcf_id}', tissue=config['tissues']),
    output:
        df_all = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['all'])
        # df_all = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai_splicemap_ref_psi']['all']
    run:
        import pandas as pd
        from tqdm import tqdm
        df = pd.concat([pd.read_parquet(i) for i in tqdm(input.df_tissue)])
        df.to_parquet(output.df_all, index=False, partition_cols=['tissue'])

# CAT infer (single CAT)     
if 'cat_infer_single' in config['models']:       
    rule splicing_pred_infer_cat:
        input:
            var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
            mmsplice_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap'],
            splicemap_5 = expand(config_static['junction_annotation']['splicemap']['psi5'],
                                genome=config['genome'], tissue=config['tissues']),
            splicemap_3 = expand(config_static['junction_annotation']['splicemap']['psi3'],
                                genome=config['genome'], tissue=config['tissues']),
            cat_count_table = expand(OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['count_table']['updated'],
                                    tissue='{tissue_cat}', genome=config['genome'])[0],
        params:
            tissue_cat = '{tissue_cat}',
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 64000,
            threads = 1
        output:
            result = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['single_cat'])
            # result = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['single_cat']
        script:
            "../infer_cat.py"
        
# CAT infer (all CATs)
if 'cat_infer_all' in config['models']:
    rule splicing_pred_mmsplice_splicemap_cat_combine:
        input:
            pred_mmsplice_cat = expand(OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['single_cat'],
                                    tissue_cat=config['tissues_cat'], vcf_id='{vcf_id}'),
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 64000,
            threads = 1
        output:
            result = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['all_cats'])
            # result = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['all_cats']
        run:
            df_mmsplice_cat = pd.concat([pd.read_parquet(i).reset_index() for i in input.pred_mmsplice_cat])
            df_mmsplice_cat.to_parquet(output.result, partition_cols='tissue_cat', index=False)
        
        
# --------------------------------AbSplice input------------------------------- 
# AbSplice-DNA (input)
if 'absplice_dna' in config['models'] or 'absplice_dna_input' in config['models']:
    rule splicing_pred_absplice_dna_input:
        input:
            pred_mmsplice = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap'],
            pred_spliceai = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'],
            gene_map = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
            gene_tpm = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_tpm'],
            var_samples_df = OUTPUT_DIR_VCF_ANNO + config_static['vcf_annotation']['rare_vars']['qual_filtered'],
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 64000,
            threads = 1
        output:
            absplice_dna_input = directory(ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna'])
            # absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna']
        script:
            "../absplice_dna_input.py"
        
# AbSplice-RNA single CAT (input)    
if 'absplice_rna_single' in config['models'] or 'absplice_rna_single_input' in config['models']:
    rule splicing_pred_absplice_rna_single_cat_input:
        input:
            absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna'],
            pred_mmsplice_cat = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['single_cat'],
            CAT_pval = expand(OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['tissue_cat_pval'],
                            tissue='{tissue_cat}', vcf_id='{vcf_id}', 
                            fraser_version='{fraser_version}', delta_psi_cutoff='{delta_psi_cutoff}', padjust_cutoff='{padjust_cutoff}', totalCounts_cutoff='{totalCounts_cutoff}'),
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 64000,
            threads = 1
        output:
            absplice_rna_input = directory(ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['rna']['single_cat'])
            # absplice_rna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['rna']['single_cat']
        script:
            "../absplice_rna_input.py"
         
# AbSplice-RNA all CATs (input)          
if 'absplice_rna_all' in config['models'] or 'absplice_rna_all_input' in config['models']:
    rule splicing_pred_absplice_rna_all_cats_input:
        input:
            absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna'],
            pred_mmsplice_cat = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap_cat']['all_cats'],
            CAT_pval = expand(OUTPUT_DIR_OUTLIER + config_static['outlier_ground_truth']['combine_gene_junction']['variant_nearest_outlier']['tissue_cat_pval'],
                            tissue=config['tissues_cat'], vcf_id='{vcf_id}',
                            fraser_version='{fraser_version}', delta_psi_cutoff='{delta_psi_cutoff}', padjust_cutoff='{padjust_cutoff}', totalCounts_cutoff='{totalCounts_cutoff}'),
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 64000,
            threads = 1
        output:
            absplice_rna_input = directory(ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['rna']['all_cats'])
            # absplice_rna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['rna']['all_cats']
        script:
            "../absplice_rna_input.py"


# --------------------------------predict AbSplice------------------------------- 
if 'absplice_dna' in config['models']:
    rule splicing_pred_absplice_dna:
        input:
            absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna'],
            absplice_model = config_precomputed['absplice']['dna'],
        params:
            median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
            tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
            extra_info = config['filtering_params']['absplice']['extra_info'],
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 32000,
            threads = 1
        output:
            absplice_dna_pred = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna'])
            # absplice_dna_pred = OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna']
        script:
            "../absplice_dna.py"


if 'pangolin_splicemaps' in config['models']:
    rule splicing_pred_pangolin_splicemaps:
        input:
            pred_pangolin = OUTPUT_DIR_SPLICING + config_static['splicing_pred']['levels']['postprocess'] + config_static['splicing_pred']['models']['pangolin'],
            splicemap_5 = expand(config_static['junction_annotation']['splicemap']['psi5'],
                             genome=config['genome'], tissue=config['tissues']),
            splicemap_3 = expand(config_static['junction_annotation']['splicemap']['psi3'],
                                genome=config['genome'], tissue=config['tissues']),
        params:
            subset_cols = config['subset_cols_absplice_models_new']
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 128000,
            threads = 1
        output:
            pangolin_splicemaps = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['pangolin_splicemaps'])
        script:
            "../pangolin_join_splicemap_info.py"


if 'absplice_dna_custom_new' in config['models']:
    rule splicing_pred_absplice_dna_custom_new:
        input:
            pangolin_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['pangolin_splicemaps'],
            mmsplice_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap'],
            splicemap_5 = expand(config_static['junction_annotation']['splicemap']['psi5'],
                             genome=config['genome'], tissue=config['tissues']),
            splicemap_3 = expand(config_static['junction_annotation']['splicemap']['psi3'],
                                genome=config['genome'], tissue=config['tissues']),
        params:
            absplice_models=config_static['absplice_models_new'],
            subset_cols = config['subset_cols_absplice_models_new']
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 128000,
            threads = 1
        output:
            absplice_dna_pred = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['absplice_custom_new']['dna'])
        script:
            "../absplice2_pangolin_tissue_mmsplice_splicemap.py"


if 'absplice2_input' in config['models']:
    rule splicing_pred_absplice2_input:
        input:
            pangolin_splicemaps = OUTPUT_DIR + config_static['splicing_pred']['models']['pangolin_splicemaps'],
            mmsplice_splicemap = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap'],
            splicemap_5 = expand(config_static['junction_annotation']['splicemap']['psi5'],
                             genome=config['genome'], tissue=config['tissues']),
            splicemap_3 = expand(config_static['junction_annotation']['splicemap']['psi3'],
                                genome=config['genome'], tissue=config['tissues']),
        params:
            subset_cols = config['subset_cols_absplice_models_new']
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 128000,
            threads = 1
        output:
            absplice2_input = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['absplice2_input'])
        script:
            "../absplice2_input.py"
        
        
if 'absplice_rna_single' in config['models']:
    rule splicing_pred_absplice_rna_single_cat:
        input:
            absplice_rna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['rna']['single_cat'],
            absplice_model = config_precomputed['absplice']['rna'],
        params:
            median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
            tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 32000,
            threads = 1
        output:
            absplice_rna_pred = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['rna']['single_cat'])
            # absplice_rna_pred = OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['rna']['single_cat']
        script:
            "../absplice_rna.py"
        

if 'absplice_rna_all' in config['models']:        
    rule splicing_pred_absplice_rna_all_cats:
        input:
            absplice_rna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['rna']['all_cats'],
            absplice_model = config_precomputed['absplice']['rna'],
        params:
            median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
            tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 32000,
            threads = 1
        output:
            absplice_rna_pred = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['rna']['all_cats'])
            # absplice_rna_pred = OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['rna']['all_cats']
        script:
            "../absplice_rna.py"


# --------------------------------AbSplice-DNA (no samples))-------------------------------     
if config['use_gene_tpm'] == True:
    gene_tpm_used = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_tpm'],
else:
    gene_tpm_used = []

if 'absplice_dna_no_samples' in config['models']: 
    rule splicing_pred_absplice_dna_input_no_samples:
        input:
            pred_mmsplice = OUTPUT_DIR + config_static['splicing_pred']['models']['mmsplice_splicemap'],
            pred_spliceai = OUTPUT_DIR + config_static['splicing_pred']['models']['spliceai'],
            gene_map = OUTPUT_DIR_JUNCTION_ANNO + config_static['junction_annotation']['gene_map'],
            gene_tpm = gene_tpm_used,
            splicemap_5 = expand(config_static['junction_annotation']['splicemap']['psi5'],
                                genome=config['genome'], tissue=config['tissues']),
            splicemap_3 = expand(config_static['junction_annotation']['splicemap']['psi3'],
                                genome=config['genome'], tissue=config['tissues']),
        params:
            tissues = config['tissues']
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 48000,
            threads = 1
        output:
            absplice_dna_input = directory(ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna_no_samples'])
            # absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna_no_samples']
        script:
            "../absplice_dna_input_no_samples.py"
        

if 'absplice_dna_no_samples' in config['models']: 
    rule splicing_pred_absplice_dna_no_samples:
        input:
            absplice_dna_input = ABSPLICE_INPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna_no_samples'],
            absplice_model = config_precomputed['absplice']['dna'],
        params:
            median_n_cutoff = config['filtering_params']['absplice']['median_n_cutoff'],
            tpm_cutoff = config['filtering_params']['absplice']['tpm_cutoff'],
            extra_info = config['filtering_params']['absplice']['extra_info'],
        resources:
            mem_mb = lambda wildcards, attempt: attempt * 48000,
            threads = 1
        output:
            absplice_dna_pred = directory(OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna_no_samples'])
            # absplice_dna_pred = OUTPUT_DIR + config_static['splicing_pred']['models']['absplice']['dna_no_samples']
        script:
            "../absplice_dna.py"


list_outputs = list()
if 'mmsplice' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_mmsplice.output, 
        vcf_id=wildcard_vcf_id), 
    )
if 'mtsplice' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_mtsplice.output, 
        vcf_id=wildcard_vcf_id), 
    )
if 'spliceai' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_spliceai.output, 
        vcf_id=wildcard_vcf_id), 
    )
if 'pangolin' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_pangolin.output, 
        vcf_id=wildcard_vcf_id), 
    )
if 'cadd_splice' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_cadd_splice.output, 
        vcf_id=wildcard_vcf_id), 
    )
if 'mmsplice_splicemap' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_mmsplice_splicemap.output, 
        vcf_id=wildcard_vcf_id), 
    )
if 'mmsplice_splicemap_ref_psi' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_mmsplice_splicemap.output, 
        vcf_id=wildcard_vcf_id),
    )
if 'spliceai_splicemap' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_spliceai_splicemap_all.output, 
        vcf_id=wildcard_vcf_id),
    )
if 'spliceai_splicemap_ref_psi' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_spliceai_splicemap_ref_psi_all.output, 
        vcf_id=wildcard_vcf_id),
    )
if 'cat_infer_single' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_infer_cat.output, 
        vcf_id=wildcard_vcf_id, tissue_cat=config['tissues_cat']),
    )
if 'cat_infer_all' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_mmsplice_splicemap_cat_combine.output, 
        vcf_id=wildcard_vcf_id),
    )
if 'absplice_dna__input' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_dna_input.output, 
        vcf_id=wildcard_vcf_id),
    )
if 'absplice_dna_no_samples__input' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_dna_input_no_samples.output, 
        vcf_id=wildcard_vcf_id),
    )
if 'absplice_rna_single__input' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_rna_single_cat_input.output, 
        vcf_id=wildcard_vcf_id, tissue_cat=config['tissues_cat']),
    )
if 'absplice_rna_all__input' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_rna_all_cats_input.output, 
        vcf_id=wildcard_vcf_id),
    )
if 'absplice_dna' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_dna.output, 
        vcf_id=wildcard_vcf_id),
    )
if 'absplice_dna_custom_new' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_dna_custom_new.output, 
        vcf_id=wildcard_vcf_id),
    )
if 'absplice_dna_custom_old' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_dna_custom_old.output, 
        vcf_id=wildcard_vcf_id),
    )
if 'absplice_dna_custom_new_no_samples' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_dna_custom_new_no_samples.output, 
        vcf_id=wildcard_vcf_id),
    )
if 'pangolin_splicemaps' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_pangolin_splicemaps.output, 
        vcf_id=wildcard_vcf_id),
    )
if 'absplice2_input' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice2_input.output, 
        vcf_id=wildcard_vcf_id),
    )
if 'absplice_dna_custom_old_no_samples' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_dna_custom_old_no_samples.output, 
        vcf_id=wildcard_vcf_id),
    )
if 'absplice_dna_no_samples' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_dna_no_samples.output, 
        vcf_id=wildcard_vcf_id),
    )
if 'absplice_rna_single' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_rna_single_cat.output, 
        vcf_id=wildcard_vcf_id, tissue_cat=config['tissues_cat']),
    )
if 'absplice_rna_all' in config['models']:
    list_outputs.append(
        expand(rules.splicing_pred_absplice_rna_all_cats.output, 
        vcf_id=wildcard_vcf_id),
    )


rule all_splicing_pred:
    input:
        list_outputs
        
        
del OUTPUT_DIR
del ABSPLICE_INPUT_DIR