# Frequently Asked Questions

## How do I run ARETE in a Slurm HPC environment?

- Set a config file under `~/.nextflow/config` to use the slurm executor:

```json
  process {
    executor = 'slurm'
    pollInterval = '60 sec'
    submitRateLimit = '60/1min'
    queueSize = 100
    // If an account is necessary:
    clusterOptions = '--account=<my-account>'
  }
```

See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#scope-executor) for a description of these options.

- Now, when running ARETE, you'll need to set additional options if your compute nodes don't have network access - as is common for most Slurm clusters. The example below uses the default test data, i.e. the `test` profile, for demonstration purposes only.

```bash
nextflow run beiko-lab/ARETE \
  --db_cache path/to/db_cache \
  --bakta_db path/to/baktadb \
  -profile test,singularity
```

Apart from `-profile singularity`, which just makes ARETE use Singularity/Apptainer containers for running the tools, there are two additional parameters:

- `--db_cache` should be the location for the pre-downloaded databases used in the DIAMOND alignments (i.e. Bacmet, VFDB, ICEberg2 and CAZy FASTA files) and in the Kraken2 taxonomic read classification.

      - Although these tools run by default, you can change the selection of annotation tools by changing `--annotation_tools` and
        skip Kraken2 by adding `--skip_kraken`. See the [parameter documentation](https://beiko-lab.github.io/arete/params/) for a full list of parameters and their defaults.

- `--bakta_db` should be the location of the pre-downloaded [Bakta database](https://github.com/oschwengers/bakta#database-download)

      - Alternatively, you can use Prokka for annotating your assemblies, since it doesn't require a downloaded database (`--use_prokka`).
