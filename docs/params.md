# beiko-lab/ARETE pipeline parameters

AMR/VF LGT-focused bacterial genomics workflow

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input_sample_table` | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row.</small></details>| `string` |  | True |  |
| `outdir` | Path to the output directory where the results will be saved. | `string` | ./results |  |  |
| `db_cache` | Directory where the databases are located | `string` |  |  |  |
| `email` | Email address for completion summary. <details><summary>Help</summary><small>Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.</small></details>| `string` |  |  |  |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |  |  |

## Reference genome options

Reference and outgroup genome fasta files required for the workflow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `reference_genome` | Path to FASTA reference genome file. | `string` |  |  |  |

## QC



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `run_checkm` | Run CheckM QC software | `boolean` |  |  |  |
| `apply_filtering` | Filter assemblies on QC results | `boolean` |  |  |  |
| `skip_kraken` | Don't run Kraken2 taxonomic classification | `boolean` |  |  |  |
| `min_n50` | Minimum N50 for filtering | `integer` | 10000 |  |  |
| `min_contigs_1000_bp` | Minimum number of contigs with >1000bp | `integer` | 1 |  |  |
| `min_contig_length` | Minimum average contig length | `integer` | 1 |  |  |

## Annotation

Parameters for the annotation subworkflow

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `annotation_tools` | Comma-separated list of annotation tools to run | `string` | mobsuite,rgi,cazy,vfdb,iceberg,bacmet,islandpath,phispy,report |  |  |
| `bakta_db` | Path to the BAKTA database | `string` |  |  |  |
| `use_prokka` | Use Prokka (not Bakta) for annotating assemblies | `boolean` |  |  |  |
| `min_pident` | Minimum match identity percentage for filtering | `integer` | 60 |  |  |
| `min_qcover` | Minimum coverage of each match for filtering | `number` | 0.6 |  |  |
| `skip_profile_creation` | Skip annotation feature profile creation | `boolean` |  |  |  |
| `feature_profile_columns` | Columns to include in the feature profile | `string` | mobsuite,rgi,cazy,vfdb,iceberg,bacmet |  |  |

## Phylogenomics

Parameters for the phylogenomics subworkflow

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `skip_phylo` | Skip Pangenomics and Phylogenomics subworkflow | `boolean` |  |  |  |
| `use_ppanggolin` | Use ppanggolin for calculating the pangenome | `boolean` |  |  |  |
| `use_full_alignment` | Use full alignment | `boolean` |  |  |  |
| `use_fasttree` | Use FastTree | `boolean` | True |  |  |

## PopPUNK

Parameters for the lineage subworkflow

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `skip_poppunk` | Skip PopPunk | `boolean` |  |  |  |
| `poppunk_model` | Which PopPunk model to use (bgmm, dbscan, refine, threshold or lineage) | `string` |  |  |  |
| `run_poppunk_qc` | Whether to run the QC step for PopPunk | `boolean` |  |  |  |
| `enable_subsetting` | Enable subsetting workflow based on genome similarity | `boolean` |  |  |  |
| `core_similarity` | Similarity threshold for core genomes | `number` | 99.99 |  |  |
| `accessory_similarity` | Similarity threshold for accessory genes | `number` | 99 |  |  |

## Recombination

Parameters for the recombination subworkflow

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `run_recombination` | Run Recombination | `boolean` |  |  |  |
| `run_verticall` | Run Verticall recombination tool | `boolean` | True |  |  |
| `run_gubbins` | Run Gubbins recombination tool | `boolean` |  |  |  |

## Dynamics



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `run_evolccm` | Run the community coevolution model | `boolean` |  |  |  |
| `run_rspr` | Run rSPR | `boolean` |  |  |  |
| `min_rspr_distance` | Minimum rSPR distance used to define processing groups | `integer` | 10 |  |  |
| `min_branch_length` | Minimum rSPR branch length | `integer` | 0 |  |  |
| `max_support_threshold` | Maximum rSPR support threshold | `integer` | 0 |  |  |
| `max_approx_rspr` | Maximum approximate rSPR distance for filtering | `integer` | -1 |  |  |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` | https://raw.githubusercontent.com/nf-core/configs/master |  | True |
| `hostnames` | Institutional configs hostname. | `string` |  |  | True |
| `config_profile_name` | Institutional config name. | `string` |  |  | True |
| `config_profile_description` | Institutional config description. | `string` |  |  | True |
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |

## Max job request options

Set the top limit for requested resources for any single job.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `max_cpus` | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`</small></details>| `integer` | 16 |  | True |
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`</small></details>| `string` | 128.GB |  | True |
| `max_time` | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`</small></details>| `string` | 240.h |  | True |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `help` | Display help text. | `boolean` |  |  | True |
| `publish_dir_mode` | Method used to save pipeline results to output directory. <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | True |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. <details><summary>Help</summary><small>An email address to send a summary email to when the pipeline is completed - ONLY sent if the pipeline does not exit successfully.</small></details>| `string` |  |  | True |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |  | True |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |  | True |
| `tracedir` | Directory to keep pipeline Nextflow logs and reports. | `string` | ${params.outdir}/pipeline_info |  | True |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |
| `show_hidden_params` | Show all params when using `--help` <details><summary>Help</summary><small>By default, parameters set as _hidden_ in the schema are not shown on the command line when a user runs with `--help`. Specifying this option will tell the pipeline to show all parameters.</small></details>| `boolean` |  |  | True |
| `enable_conda` | Run this workflow with Conda. You can also use '-profile conda' instead of providing this parameter. | `boolean` |  |  | True |
| `singularity_pull_docker_container` | Instead of directly downloading Singularity images for use with Singularity, force the workflow to pull and convert Docker containers instead. <details><summary>Help</summary><small>This may be useful for example if you are unable to directly pull Singularity containers to run the pipeline due to http/https proxy issues.</small></details>| `boolean` |  |  | True |
| `schema_ignore_params` |  | `string` | genomes,modules |  |  |
| `multiqc_logo` |  | `string` |  |  | True |
