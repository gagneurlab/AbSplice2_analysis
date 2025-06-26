# AbSplice2 analysis

Repository for the analyses done in the preprint (TODO: link will be added) "Aberrant splicing prediction during human development".

# Folder Structure
The project is setup as a snakemake pipeline. 

* The `workflow` folder contains the config files and the analysis scripts for each dataset.
    * `config` contains yaml files for the GTEx, mito and ALS datasets, as well as a yaml that defines the file structure used in all analysis. 
    * `scripts` contains all processing and analysis scripts. The `common` folder contains all scripts shared across datasets. `gtex_v8`, `gvex`, `devAbSplice` and `all_variants` contain analysis scripts specific to the respective analysis. Each of these folders contains a Snakefile that runs the analysis of this respective analysis. The analysis should be executed in the following order, as there are intermediate results required as inputs (e.g. AbSplice2 model, SpliceMaps from GTEx tissues/ Developmental SpliceMaps): 1) `gtex_v8`, 2) `gvex`, 3) `devAbSplice`, 4) `all_variants`

# Overview of the analysis pipeline
In `scripts/common` most of the scripts of the pipeline are included.
The main tasks of the pipeline are:
* `vcf_annotation` filters the provided vcf file for rare, high quality variants 
* `junction_annotation` generates SpliceMaps based on RNA-seq split read counts from FRASER and FRASER2
* `outlier_ground_truth` generates outliers (FRASER and FRASER2) on gene and junction level, filtered for significance and rare variants in the vicinity
* `splicing_pred` generates splicing predictions for SpliceAI, Pangolin, MMSplice, AbSplice and AbSplice2
* `splicing_result` provides maximum predictions on a gene and variant level
* `benchmark` creates the benchmark dataframes from predictions and outliers. Where each gene, sample, tissue combination in the benchmark contains at least one rare variant on the gene

# Necessary files to provide by the user
In order for the pipeline to work users need to provide some files:

* create the folder `./workflow/data/` (recommended to use a symlink pointing to sufficient storage)

* outlier calls from DROP. Create a symlink to the DROP results in the folder `./workflow/data/resources/{dataset}/DROP/`, where dataset can be `gtex_v8`, `gvex` and `devAbSplice`.  
We used DROP v.1.4.0. For GTEx v8 the necessary config file to run DROP is provided in [here](https://github.com/gagneurlab/AbSplice2_analysis/blob/master/workflow/_data/resources/gtex_v8/DROP_config.yaml): `./workflow/_data/resources/gtex_v8/DROP_config.yaml`, with the corresponding sample annotation for DROP in [here](https://github.com/gagneurlab/AbSplice2_analysis/blob/master/workflow/_data/resources/gtex_v8/DROP_sample_annotation.tsv): `./workflow/_data/resources/gtex_v8/DROP_sample_annotation.tsv`.

* vcf files from the dataset. Users need to normalize the vcf files and store them in `./workflow/data/resources/{dataset}/vcf_normalized/{vcf_id}.vcf.gz`, where dataset can be `gtex_v8` and `gvex`. `vcf_id` can be anything (e.g. sampleID or split huge vcf by chromosome to parallelize computations).

* gtf and fasta files. For each analysis create a symlink to the gtf (`./workflow/data/resources/{dataset}/gtf_file`) and fasta file (`./workflow/data/resources/{dataset}/fasta_file`) used in the respective analysis:
    * gtex_v8: gtf_file (gencode.v38/gencode.v38.annotation.gtf), fasta_file (gencode.v38/GRCh38.primary_assembly.genome.fa)
    * gvex: gtf_file (release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz), fasta_file (release_34/GRCh37_mapping/GRCh37.primary_assembly.genome.fa)
    * devAbSplice: gtf_file (release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz), fasta_file (release_34/GRCh37_mapping/GRCh37.primary_assembly.genome.fa)

* databases required during pipeline. create the folders `.workflow/data/resources/common/{genome}/`, where genome is `hg19` and `hg38`. To those folders download:
    * pangolin database. Create a symlin in `.workflow/data/resources/common/{genome}/pangolin_db`, where genome can be `hg19` and `hg38`. The database can be created from: https://github.com/tkzeng/Pangolin/tree/main
    * rocksdb database with precomputed SpliceAI scores. Create a symlink in `./workflow/data/resources/common/spliceai_db/spliceAI_{genome}_{chrom}.db`, where genome can be `hg19` and `hg38`. These databases can be created from: https://github.com/gagneurlab/spliceai_rocksdb.
    * rocksdb database with gnomAD minor allele frequencies. Create a symlink in `./workflow/data/resources/common/{genome}/gnomAD_maf_db/rocksdb/maf.db`, where genome can be `hg19` and `hg38`. These databases can be created from: https://github.com/gagneurlab/gnomad_rocksdb.

* copy the files that are provided in this github repo in the folder `./workflow/_data/` to the data folder you created (replacing `_data` with `data`)


# Run the pipeline
Follow these steps to run the pipeline:
- Install the absplice environment. Follow the steps indicated in the [absplice GitHub](https://github.com/gagneurlab/absplice/tree/master). You can use the docker image or conda environment.
- Activate the absplice environment.
- From the absplice environment, run the snakemake workflow for `gtex_v8` (below just with 1 job; increase this depending on your available resources):
    ```
    cd workflow/scripts/gtex_v8
    python -m snakemake -j 1 --use-conda
    ```
- Similarly run the snakemake workflow for `gvex`, `devAbSplice` and `all_variants`. 

# Figure generation
All figures from the manuscript can be generated from the the provided scripts in [here](https://github.com/gagneurlab/AbSplice2_analysis/tree/master/figures_R). The minimal dataset to produce the figures is available (TODO: will be provided). To run the scripts, first download and unzip the provided data folder to the root directory of this repository.








