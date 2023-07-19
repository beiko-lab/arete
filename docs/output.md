# beiko-lab/ARETE: Output

## Introduction

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps (steps in _italics_ don't run by default):

- [FastQC](#fastqc) - Raw and trimmed read QC
- [FastP](#fastp) - Read trimming
- [Kraken2](#kraken2) - Taxonomic assignment
- [Unicycler](#unicycler) - Short read assembly
- [_Quast_](#quast) - Assembly quality score
- [Bakta](#bakta) or [_Prokka_](#prokka) - Gene detection and annotation
- [Panaroo](#panaroo) - Pangenome alignment
- [MobRecon](#mobrecon) - Reconstruction and typing of plasmids
- [RGI](#rgi) - Detection and annotation of AMR determinants
- [IslandPath](#islandpath) - Predicts genomic islands in bacterial and archaeal genomes.
- [PhiSpy](#phispy) - Prediction of prophages from bacterial genomes
- [_IntegronFinder_](#integronfinder) - Finds integrons in DNA sequences
- [Diamond](#diamond) - Detection and annotation of genes using external databases.
  - CAZy: Carbohydrate metabolism
  - VFDB: Virulence factors
  - BacMet: Metal resistance determinants
  - ICEberg: Integrative and conjugative elements
- [PopPUNK](#poppunk) - Genome clustering
- [_Verticall_](#verticall) - Conduct pairwise assembly comparisons between genomes in a same PopPUNK cluster
- [_SKA2_](#ska2) - Generate a whole-genome FASTA alignment for each genome within a cluster.
- [_Gubbins_](#gubbins) - Detection of recombination events within genomes of the same cluster.
- [FastTree](#fasttree) or [_IQTree_](#iqtree) - Maximum likelihood core genome phylogenetic tree
- [_SNPsites_](#snpsites) - Extracts SNPs from a multi-FASTA alignment
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

<!-- TODO put all the other crap in below. Can't be arsed today -->

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `read_processing/*_fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics for your untrimmed raw fastq files.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

**NB:** The FastQC plots in this directory are generated relative to the raw, input reads. They may contain adapter sequence and regions of low quality. To see how your reads look after adapter and quality trimming please refer to the FastQC reports in the `trimgalore/fastqc/` directory.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### fastp

<details markdown="1">
<summary>Output files</summary>

- `read_processing/fastp/`
  - `${meta.id}` : Trimmed files and trimming reports for each input sample.

</details>

[fastp](https://github.com/OpenGene/fastp) is a all-in-one fastq preprocessor for read/adapter trimming and quality control. It is used in this pipeline for trimming adapter sequences and discard low-quality reads.

### Kraken2

<details markdown="1">
<summary>Output files</summary>

- `read_processing/kraken2/`
  - `*.kraken2.report.txt` : Text file containing genome-wise information of Kraken2 findings. See [here](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#output-formats) for details.
  - `*.classified(_(1|2))?.fastq.gz` : Fasta file containing classified reads. If paired-end, one file per end.
  - `*.unclassified(_(1|2))?.fastq.gz` : Fasta file containing unclassified reads. If paired-end, one file per end.

</details>

Kraken2 is a read classification software which will assign taxonomy to each read comprising a sample. These results may be analyzed as an indicator of contamination.

### Unicycler

<details markdown="1">
<summary>Output files</summary>

- `assembly/unicycler/`
  - `*.assembly.gfa`
  - `*.scaffolds.fa`
  - `*.unicycler.log`

</details>

Short/hybrid read assembler. For now only handles short reads in ARETE.

### Quast

<details markdown="1">
<summary>Output files</summary>

- `assembly/quast/`
  - `report.tsv` : A tab-seperated report compiling all QC metrics recorded over all genomes
  - `quast/`
    - `report.(html|tex|pdf|tsv|txt)`: The Quast report in different file formats
    - `transposed_report.(tsv|txt)` : Transpose of the Quast report
    - `quast.log` : Log file of all Quast runs
    - `icarus_viewers/`
      - `contig_size_viewer.html`
    - `basic_stats/`: Directory containing various summary plots generated by Quast.

</details>

### Bakta

<details markdown="1">
<summary>Output files</summary>

- `annotation/bakta/`
  - `${sample_id}/` : Bakta results will be in one directory per genome.
    - `${sample_id}.tsv`: annotations as simple human readble TSV
    - `${sample_id}.gff3`: annotations & sequences in GFF3 format
    - `${sample_id}.gbff`: annotations & sequences in (multi) GenBank format
    - `${sample_id}.embl`: annotations & sequences in (multi) EMBL format
    - `${sample_id}.fna`: replicon/contig DNA sequences as FASTA
    - `${sample_id}.ffn`: feature nucleotide sequences as FASTA
    - `${sample_id}.faa`: CDS/sORF amino acid sequences as FASTA
    - `${sample_id}.hypotheticals.tsv`: further information on hypothetical protein CDS as simple human readble tab separated values
    - `${sample_id}.hypotheticals.faa`: hypothetical protein CDS amino acid sequences as FASTA
    - `${sample_id}.txt`: summary as TXT
    - `${sample_id}.png`: circular genome annotation plot as PNG
    - `${sample_id}.svg`: circular genome annotation plot as SVG

</details>

Bakta is a tool for the rapid & standardized annotation of bacterial genomes and plasmids from both isolates and MAGs

### _Prokka_

<details markdown="1">
<summary>Output files</summary>

- `annotation/prokka/`
  - `${sample_id}/` : Prokka results will be in one directory per genome.
    - `${sample_id}.err` : Unacceptable annotations
    - `${sample_id}.faa` : Protein FASTA file of translated CDS sequences
    - `${sample_id}.ffn` : Nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA)
    - `${sample_id}.fna` : Nucleotide FASTA file of input contig sequences
    - `${sample_id}.fsa` : Nucleotide FASTA file of the input contig sequences, used by "tbl2asn" to create the .sqn file. It is mostly the same as the .fna file, but with extra Sequin tags in the sequence description lines.
    - `${sample_id}.gff` : This is the master annotation in GFF3 format, containing both sequences and annotations.
    - `${sample_id}.gbk` : This is a standard Genbank file derived from the master .gff.
    - `${sample_id}.log` : Contains all the output that Prokka produced during its run. This is a record of what settings used, even if the --quiet option was enabled.
    - `${sample_id}.sqn` : An ASN1 format "Sequin" file for submission to Genbank. It needs to be edited to set the correct taxonomy, authors, related publication etc.
    - `${sample_id}.tbl` : Feature Table file, used by "tbl2asn" to create the .sqn file.
    - `${sample_id}.tsv` : Tab-separated file of all features: locus_tag,ftype,len_bp,gene,EC_number,COG,product
    - `${sample_id}.txt` : Statistics relating to the annotated features found.

</details>

Prokka is a software tool to annotate bacterial, archaeal and viral genomes quickly and produce standards-compliant output files.

### RGI

<details markdown="1">
<summary>Output files</summary>

- `annotation/rgi/`
  - `${sample_id}_rgi.txt` : A TSV report containing all AMR predictions for a given genome. For more info see [here](https://github.com/arpcard/rgi#rgi-main-tab-delimited-output-details)

</details>

RGI predicts AMR determinants using the CARD ontology and various trained models.

### MobRecon

<details markdown="1">
<summary>Output files</summary>

- `annotation/mob_recon`
  - `${sample_id}_mob_recon/` : MobRecon results will be in one directory per genome.
    - `contig_report.txt` - This file describes the assignment of the contig to chromosome or a particular plasmid grouping.
    - `mge.report.txt` - Blast HSP of detected MGE's/repetitive elements with contextual information.
    - `chromosome.fasta` - Fasta file of all contigs found to belong to the chromosome.
    - `plasmid_*.fasta` - Each plasmid group is written to an individual fasta file which contains the assigned contigs.
    - `mobtyper_results` - Aggregate MOB-typer report files for all identified plasmid.

</details>

MobRecon reconstructs individual plasmid sequences from draft genome assemblies using the clustered plasmid reference databases

### DIAMOND

<details markdown="1">
<summary>Output files</summary>

- `annotation/(vfdb|bacmet|cazy|iceberg2)/`
  - `${sample_id}/${sample_id}_(VFDB|BACMET|CAZYDB|ICEberg2).txt` : Blast6 formatted TSVs indicating BlastX results of the genes from each genome against VFDB, BacMet, and CAZy databases.
  - `(VFDB|BACMET|CAZYDB|ICEberg2).txt` : Table with all hits to this database, with a column describing which genome the match originates from. Sorted and filtered by the match's coverage.

</details>

[Diamond](https://github.com/bbuchfink/diamond) is a sequence aligner for protein and translated DNA searches, designed for high performance analysis of big sequence data. We use DIAMOND to predict the presence of virulence factors, heavy metal resistance determinants, carbohydrate-active enzymes, and integrative and conjugative elements using [VFDB](http://www.mgc.ac.cn/VFs/), [BacMet](http://bacmet.biomedicine.gu.se/), [CAZy](http://www.cazy.org/), and [ICEberg2](https://bioinfo-mml.sjtu.edu.cn/ICEberg2/index.php) respectively.

### IslandPath

<details markdown="1">
<summary>Output files</summary>

- `annotation/islandpath/`
  - `${sample_id}/` : IslandPath results will be in one directory per genome.
    - `${sample_id}.tsv` : IslandPath results
    - `Dimob.log` : IslandPath execution log

</details>

[IslandPath](https://github.com/brinkmanlab/islandpath) is a standalone software to predict genomic islands
in bacterial and archaeal genomes based on the presence of dinucleotide biases and mobility genes.

### IntegronFinder

Disabled by default. Enable by adding `--run_integronfinder` to your command.

<details markdown="1">
<summary>Output files</summary>

- `annotation/integron_finder/`
  - `Results_Integron_Finder_${sample_id}/` : IntegronFinder results will be in one directory per genome.

</details>

[Integron Finder](https://github.com/gem-pasteur/Integron_Finder) is a bioinformatics tool to find integrons in bacterial genomes.

### PhiSpy

<details markdown="1">
<summary>Output files</summary>

- `annotation/phispy/`
  - `${sample_id}/` : PhiSpy results will be in one directory per genome.

See the [PhiSpy documentation](https://github.com/linsalrob/PhiSpy#output-files) for an extensive description of the output.

</details>

[PhiSpy](https://github.com/linsalrob/PhiSpy) is a tool for identification of prophages in Bacterial (and probably Archaeal) genomes. Given an annotated genome it will use several approaches to identify the most likely prophage regions.

### Panaroo

<details markdown="1">
<summary>Output files</summary>

- `pangenomics/panaroo/results/`

See [the panaroo documentation](https://gtonkinhill.github.io/panaroo/#/gettingstarted/output) for an extensive description of output provided.

</details>

[Panaroo](https://gtonkinhill.github.io/panaroo/) is a Bacterial Pangenome Analysis Pipeline.

### FastTree

<details markdown="1">
<summary>Output files</summary>

- `phylogenomics/fasttree/`
  - `*.tre` : Newick formatted maximum likelihood tree of core-genome alignment.

</details>

[FastTree](http://www.microbesonline.org/fasttree/) infers approximately-maximum-likelihood phylogenetic trees from alignments of nucleotide or protein sequences

### _IQTree_

<details markdown="1">
<summary>Output files</summary>

- `phylogenomics/iqtree/`
  - `*.treefile` : Newick formatted maximum likelihood tree of core-genome alignment.

</details>

[IQTree](http://www.iqtree.org/) is a fast and effective stochastic algorithm to infer phylogenetic trees by maximum likelihood.

### _SNPsites_

<details markdown="1">
<summary>Output files</summary>

- `phylogenomics/snpsites/`
  - `filtered_alignment.fas` : Variant fasta file.
  - `constant.sites.txt` : Text file containing counts of constant sites.

</details>

[SNPsites](https://github.com/sanger-pathogens/snp-sites) is a tool to rapidly extract SNPs from a multi-FASTA alignment.

### PopPUNK

<details markdown="1">
<summary>Output files</summary>

- `poppunk_results/`
  - `poppunk_${poppunk_model}/` - Results from PopPUNK's fit-model command
  - `poppunk_visualizations/` - Results from the poppunk_visualise command

</details>

[PopPUNK](https://poppunk.net/) is a tool for clustering genomes.

### Verticall

<details markdown="1">
<summary>Output files</summary>

- `recombination/verticall/`
  - `verticall_cluster*.tsv` - Verticall results for the genomes within this PopPUNK cluster.

</details>

[Verticall](https://github.com/rrwick/Verticall/) is a tool to help produce bacterial genome phylogenies which are not influenced by horizontally acquired sequences such as recombination.

### SKA2

<details markdown="1">
<summary>Output files</summary>

- `recombination/ska2/`
  - `cluster_*.aln` - SKA2 results for the genomes within this PopPUNK cluster.

</details>

[SKA2](https://github.com/bacpop/ska.rust) (Split Kmer Analysis) is a toolkit for prokaryotic (and any other small, haploid) DNA sequence analysis using split kmers.

### Gubbins

<details markdown="1">
<summary>Output files</summary>

- `recombination/gubbins/`
  - `cluster_*/` - Gubbins results for the genomes within this PopPUNK cluster.
  </details>

[Gubbins](https://github.com/nickjcroucher/gubbins) is an algorithm that iteratively identifies loci containing elevated densities of base substitutions while concurrently constructing a phylogeny based on the putative point mutations outside of these regions.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
