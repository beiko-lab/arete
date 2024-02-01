# beiko-lab/ARETE: Usage

## Introduction

The ARETE pipeline can is designed as an end-to-end workflow manager for genome assembly, annotation, and phylogenetic analysis, beginning with read data. However, in some cases a user may wish to stop the pipeline prior to annotation or use the annotation features of the work flow with pre-existing assemblies. Therefore, ARETE allows users different use cases:

1. Run the full pipeline end-to-end.
2. Input a set of reads and stop after assembly.
3. Input a set of assemblies and perform QC.
4. Input a set of assemblies and perform annotation and taxonomic analyses.
5. Input a set of assemblies and perform genome clustering with PopPUNK.
6. Input a set of assemblies and perform phylogenomic and pangenomic analysis.

This document will describe how to perform each workflow.

["Running the pipeline"](#running-the-pipeline) will show some example command on how to use these different entries to ARETE.

## Samplesheet input

No matter your use case, you will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location.
For full runs and assembly, it has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input_sample_table '[path to samplesheet file]'
```

### Full workflow or assembly samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 4 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```bash
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,
```

| Column    | Description                                                                                                                |
| --------- | -------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample.              |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |

An [example samplesheet](https://github.com/beiko-lab/arete/blob/master/test/test_dataset.csv) has been provided with the pipeline.

### Annotation only samplesheet

The ARETE pipeline allows users to provide pre-existing assemblies to make use of the annotation and reporting features of the workflow. Users may use the `assembly_qc` entry point to perform QC on the assemblies. **Note that the QC workflow does not automatically filter low quality assemblies, it simply generates QC reports!** `annotation`, `assembly_qc` and `poppunk` workflows accept the same format of sample sheet.

The sample sheet must be a 2 column, comma-seperated CSV file with header.

| Column          | Description                                                                                                   |
| --------------- | ------------------------------------------------------------------------------------------------------------- |
| `sample`        | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. |
| `fna_file_path` | Full path to fna file for assembly or genome. File must have `.fna` file extension.                           |

An [example samplesheet](https://github.com/beiko-lab/arete/blob/master/test/test_annotation_dataset.csv) has been provided with the pipeline.

### Phylogenomics and Pangenomics only samplesheet

The ARETE pipeline allows users to provide pre-existing assemblies to make use of the phylogenomic and pangenomic features of the workflow.

The sample sheet must be a 2 column, comma-seperated CSV file with header.

| Column          | Description                                                                                                                                                                                |
| --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `sample`        | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample.                                                                              |
| `gff_file_path` | Full path to GFF file for assembly or genome. File must have `.gff` or `.gff3` file extension. These files can be the ones generated by Prokka or Bakta in ARETE's annotation subworkflow. |

## Reference Genome

For full workflow or assembly, users may provide a path to a reference genome in fasta format for use in assembly evaluation.

```bash
--reference_genome ref.fasta
```

<!--
The pipeline also requires a genome in fasta format to be supplied to use as an outgroup for phylogenetic analyses:
```bash
--outgroup_genome out.fasta
```
-->

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run beiko-lab/ARETE --input_sample_table samplesheet.csv --reference_genome ref.fasta --poppunk_model bgmm  -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

As written above, the pipeline also allows users to execute only assembly or only annotation.

### Assembly Entry

To execute assembly (reference genome optional):

```bash
nextflow run beiko-lab/ARETE -entry assembly --input_sample_table samplesheet.csv --reference_genome ref.fasta  -profile docker
```

### Assembly QC Entry

To execute QC on pre-existing assemblies (reference genome optional):

```bash
nextflow run beiko-lab/ARETE -entry assembly_qc --input_sample_table samplesheet.csv --reference_genome ref.fasta  -profile docker
```

### Annotation Entry

To execute annotation of pre-existing assemblies (PopPUNK model can be either bgmm, dbscan, refine, threshold or lineage):

```bash
nextflow run beiko-lab/ARETE -entry annotation --input_sample_table samplesheet.csv --poppunk_model bgmm -profile docker
```

### PopPUNK Entry

To execute annotation of pre-existing assemblies (PopPUNK model can be either bgmm, dbscan, refine, threshold or lineage):

```bash
nextflow run beiko-lab/ARETE -entry poppunk --input_sample_table samplesheet.csv --poppunk_model bgmm -profile docker
```

### Phylogenomics and Pangenomics Entry

To execute phylogenomic and pangenomics analysis on pre-existing assemblies:

```bash
nextflow run beiko-lab/ARETE -entry phylogenomics --input_sample_table samplesheet.csv -profile docker
```

### rSPR Entry

To execute the rSPR analysis on pre-existing trees:

```bash
nextflow run beiko-lab/ARETE \
  -entry rspr \
  --input_sample_table samplesheet.csv \
  --core_gene_tree core_gene_alignment.tre \
  --concatenated_annotation BAKTA.txt \
  -profile docker
```

The parameters being:

- `--core_gene_tree` - The reference tree, coming from a core genome alignment, like the one generated by panaroo in ARETE.
- `--concatenated_annotation` - The tabular annotation results (TSV) for all genomes, like the ones generated at the end of Prokka or Bakta in ARETE. Although useful, it's not necessary to execute the rSPR entry.
- `--input_sample_table` - A samplesheet containing all individual gene trees in the following format:

`gene_tree,path
CDS_0000,/path/to/CDS_0000.tre
CDS_0001,/path/to/CDS_0001.tre
CDS_0002,/path/to/CDS_0002.tre
CDS_0003,/path/to/CDS_0003.tre
CDS_0004,/path/to/CDS_0004.tre
`

### evolCCM Entry

To execute the evolCCM analysis on a pre-existing reference tree and feature profile:

```bash
nextflow run beiko-lab/ARETE \
  -entry evolccm \
  --core_gene_tree core_gene_alignment.tre \
  --feature_profile feature_profile.tsv.gz \
  -profile docker
```

The parameters being:

- `--core_gene_tree` - The reference tree, coming from a core genome alignment, like the one generated by panaroo in ARETE.
- `--feature_profile` - A presence/absence TSV matrix of features
  in genomes. Genome names should be the same in the core tree and
  should be contained to a 'genome_id' column, with all other columns represent features absent (0) or present (1) in each genome. I.e.:

```
genome_id	plasmid_AA155	plasmid_AA161
ED010	0	0
ED017	0	1
ED040	0	0
ED073	0	1
ED075	1	1
ED082	0	1
ED142	0	1
ED178	0	1
ED180	0	0
```

### Recombination Entry

To execute the recombination analysis on pre-existing assemblies (PopPUNK model can be either bgmm, dbscan, refine, threshold or lineage):

```bash
nextflow run beiko-lab/ARETE \
  -entry recombination \
  --input_sample_table samplesheet.csv \
  --poppunk_model dbscan \
  -profile docker
```

### Gene Order Entry

To execute the Gene Order analysis on pre-existing assemblies and RGI annotations:

```bash
nextflow run beiko-lab/ARETE \
  -entry gene_order \
  --input_sample_table gene_order_samplesheet.csv \
  -profile docker
```

- `--input_sample_table` - A samplesheet containing a fasta file, a genbank file and an RGI output file for each assembly:

```
sample,fna_file_path,gbk,rgi
SAMD00052607,SAMD00052607.faa,SAMD00052607.gbk,SAMD00052607_rgi.txt
SAMEA1466699,SAMEA1466699.faa,SAMEA1466699.gbk,SAMEA1466699_rgi.txt
SAMEA1486355,SAMEA1486355.faa,SAMEA1486355.gbk,SAMEA1486355_rgi.txt
```

## Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull beiko-lab/ARETE
```

## Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [ARETE releases page](https://github.com/beiko-lab/ARETE/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`

      - A generic configuration profile to be used with [Docker](https://docker.com/)

- `singularity`

      - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)

- `podman`

      - A generic configuration profile to be used with [Podman](https://podman.io/)

- `shifter`

      - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)

- `charliecloud`

      - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)

- `test`

      - A profile with a complete configuration for automated testing
      - Can run in personal computers with at least 6GB of RAM and 2 CPUs
      - Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customise the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `UNICYCLER` 32GB of memory, you could use the following config:

```nextflow
process {
  withName: UNICYCLER {
    memory = 32.GB
  }
}
```

To find the exact name of a process you wish to modify the compute resources, check the live-status of a nextflow run displayed on your terminal or check the nextflow error for a line like so: `Error executing process > 'bwa'`. In this case the name to specify in the custom config file is `bwa`.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

Sometimes LevelDB, which is used by Nextflow to track execution metadata, can lead to memory-related issues, often showing as a **SIGBUS** error. [This tends to happen when running Nextflow in SLURM environments](https://github.com/nextflow-io/nextflow/issues/842).

In this case, setting `NXF_OPTS="-Dleveldb.mmap=false"` in your `~/.bashrc` or immediately before executing `nextflow run` usually solves the issue.

## ARETE's storage requirements

ARETE generates a lot of intermediary files, which is even further exacerbated if you are running on a dataset with more than 100 genomes.
Before running ARETE you should make sure you have at least 500 GB of free storage.

After running ARETE and checking your results, you can remove the `work/` directory in your working directory, which is where Nextflow stores its cache.
**Be aware that deleting `work/` will make it so your pipeline won't re-run with cache when using the `-resume` flag, it will run every process from scratch.**
