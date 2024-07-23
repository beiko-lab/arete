# beiko-lab/ARETE pipeline parameters

AMR/VF LGT-focused bacterial genomics workflow

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `input_sample_table` | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row.</small></details>| `string` |  |
| `outdir` | Path to the output directory where the results will be saved. | `string` | ./results |
| `db_cache` | Directory where the databases are located | `string` |  |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>| `string` |  |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |

## Reference genome options

Reference and outgroup genome fasta files required for the workflow.

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `reference_genome` | Path to FASTA reference genome file. | `string` |  |

## QC



| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `run_checkm` | Run CheckM QC software | `boolean` |  |
| `apply_filtering` | Filter assemblies on QC results | `boolean` |  |
| `skip_kraken` | Don't run Kraken2 taxonomic classification | `boolean` |  |
| `min_n50` | Minimum N50 for filtering | `integer` | 10000 |
| `min_contigs_1000_bp` | Minimum number of contigs with >1000bp | `integer` | 1 |
| `min_contig_length` | Minimum average contig length | `integer` | 1 |

## Annotation

Parameters for the annotation subworkflow

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `annotation_tools` | Comma-separated list of annotation tools to run | `string` | mobsuite,rgi,cazy,vfdb,iceberg,bacmet,islandpath,phispy,report |
| `bakta_db` | Path to the BAKTA database | `string` |  |
| `use_prokka` | Use Prokka (not Bakta) for annotating assemblies | `boolean` |  |
| `min_pident` | Minimum match identity percentage for filtering | `integer` | 60 |
| `min_qcover` | Minimum coverage of each match for filtering | `number` | 0.6 |
| `skip_profile_creation` | Skip annotation feature profile creation | `boolean` |  |
| `feature_profile_columns` | Columns to include in the feature profile | `string` | mobsuite,rgi,cazy,vfdb,iceberg,bacmet |

## Phylogenomics

Parameters for the phylogenomics subworkflow

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `skip_phylo` | Skip Pangenomics and Phylogenomics subworkflow | `boolean` |  |
| `use_ppanggolin` | Use ppanggolin for calculating the pangenome | `boolean` |  |
| `use_full_alignment` | Use full alignment | `boolean` |  |
| `use_fasttree` | Use FastTree | `boolean` | True |
| `feature_dispersion_columns` | Columns from the input samplesheet to use in the feature dispersion module | `string` |  |

## PopPUNK

Parameters for the lineage subworkflow

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `skip_poppunk` | Skip PopPunk | `boolean` |  |
| `poppunk_model` | Which PopPunk model to use (bgmm, dbscan, refine, threshold or lineage) | `string` |  |
| `run_poppunk_qc` | Whether to run the QC step for PopPunk | `boolean` |  |
| `enable_subsetting` | Enable subsetting workflow based on genome similarity | `boolean` |  |
| `core_similarity` | Similarity threshold for core genomes | `number` | 99.99 |
| `accessory_similarity` | Similarity threshold for accessory genes | `number` | 99 |

## Gene Order

Parameters for the Gene Order Subworkflow

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `run_gene_order` | Whether to run the Gene Order subworkflow | `boolean` |  |
| `input_file_path` |  | `string` | /home/jvfe/dev/dalhousie/arete/test/gene-order/rgi_input.txt |
| `gene_order_percent_cutoff` | Cutoff percentage of genomes a gene should be present within to be included in extraction and subsequent analysis. Should a float between 0 and 1 (e.g., 0.25 means only genes present in a minimum of 25% of genomes are kept). | `number` | 0.25 |
| `gene_order_label_cols` | If using annotation files predicting features, list of space separated column names to be added to the gene names | `string` | None |
| `num_neighbors` | Neighborhood size to extract. Should be an even number N, such that N/2 neighbors upstream and N/2 neighbors downstream will be analyzed. | `integer` | 10 |
| `inflation` | Inflation hyperparameter value for Markov Clustering Algorithm. | `integer` | 2 |
| `epsilon` | Epsilon hyperparameter value for DBSCAN clustering. | `number` | 0.5 |
| `minpts` | Minpts hyperparameter value for DBSCAN clustering. | `integer` | 5 |
| `plot_clustering` | Create Clustering HTML Plots | `boolean` |  |

## Recombination

Parameters for the recombination subworkflow

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `run_recombination` | Run Recombination | `boolean` |  |
| `run_verticall` | Run Verticall recombination tool | `boolean` | True |
| `run_gubbins` | Run Gubbins recombination tool | `boolean` |  |

## Dynamics



| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `run_evolccm` | Run the community coevolution model | `boolean` |  |
| `run_rspr` | Run rSPR | `boolean` |  |
| `min_rspr_distance` | Minimum rSPR distance used to define processing groups | `integer` | 10 |
| `min_branch_length` | Minimum rSPR branch length | `integer` | 0 |
| `max_support_threshold` | Maximum rSPR support threshold | `number` | 0.7 |
| `max_approx_rspr` | Maximum approximate rSPR distance for filtering | `integer` | -1 |
| `min_heatmap_approx_rspr` | Minimum approximate rSPR distance used to generate heatmap | `integer` | 0 |
| `max_heatmap_approx_rspr` | Maximum approximate rSPR distance used to generate heatmap | `integer` | -1 |
| `min_heatmap_exact_rspr` | Minimum exact rSPR distance used to generate heatmap | `integer` | 0 |
| `max_heatmap_exact_rspr` | Maximum exact rSPR distance used to generate heatmap | `integer` | -1 |
| `core_gene_tree` | Core (or reference) genome tree. Used in the rSPR and evolCCM entries. | `string` |  |
| `concatenated_annotation` | TSV table of annotations for all genomes. Such as the ones generated by Bakta or Prokka in ARETE. | `string` |  |
| `feature_profile` | Feature profile TSV (A presence-absence matrix). Used in the evolCCM entry. | `string` |  |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` | https://raw.githubusercontent.com/nf-core/configs/master |
| `hostnames` | Institutional configs hostname. | `string` |  |
| `config_profile_name` | Institutional config name. | `string` |  |
| `config_profile_description` | Institutional config description. | `string` |  |
| `config_profile_contact` | Institutional config contact information. | `string` |  |
| `config_profile_url` | Institutional config URL link. | `string` |  |

## Max job request options

Set the top limit for requested resources for any single job.

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `max_cpus` | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`</small></details>| `integer` | 16 |
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`</small></details>| `string` | 128.GB |
| `max_time` | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`</small></details>| `string` | 240.h |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default |
|-----------|-----------|-----------|-----------|
| `help` | Display help text. | `boolean` |  |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>| `string` |  |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |
| `tracedir` | Directory to keep pipeline Nextflow logs and reports. | `string` | ${params.outdir}/pipeline_info |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |
| `show_hidden_params` | Show all params when using `--help` <details><summary>Help</summary><small>By default, parameters set as _hidden_ in the schema are not shown on the command line when a user runs with `--help`. Specifying this option will tell the pipeline to show all parameters.</small></details>| `boolean` |  |
| `enable_conda` | Run this workflow with Conda. You can also use '-profile conda' instead of providing this parameter. | `boolean` |  |
| `singularity_pull_docker_container` | Instead of directly downloading Singularity images for use with Singularity, force the workflow to pull and convert Docker containers instead. <details><summary>Help</summary><small>This may be useful for example if you are unable to directly pull Singularity containers to run the pipeline due to http/https proxy issues.</small></details>| `boolean` |  |
| `schema_ignore_params` |  | `string` | genomes,modules |
| `multiqc_logo` |  | `string` |  |
