## PopPUNK subsetting

The subsampling subworkflow is executed if you want
to reduce the number of genomes that get added to the phylogenomics subworkflow.
By reducing the number of genomes, you can potentially reduce resource requirements for the pangenomics and phylogenomics tools.

To enable this subworkflow, add `--enable_subsetting` when running beiko-lab/ARETE.
This will subset genomes based on their core genome similarity and accessory genome similarity, as calculated via their PopPUNK distances.

By default, the threshold is `--core_similarity 99.9` and `--accessory_similarity 99`. But these
can be changed by adding these parameters to your execution. What happens then is if any pair of genomes is this similar, only one genome from this pair will
be included in the phylogenomic section. All of the removed genome IDs will be present under `poppunk_results/removed_genomes.txt`.

By adding `--enable_subsetting`, you'll be adding two processes to the execution DAG:

- POPPUNK_EXTRACT_DISTANCES: This process will extract pair-wise distances between all genomes, returning a table under `poppunk_results/distances/`. This table will be used to perform the subsetting.
- MAKE_HEATMAP: This process will create a heatmap showing different similarity thresholds and the number of genomes that'd be present in each of the possible subsets. It'll also be under `poppunk_results/distances/`.

### Example command

The command below will execute the 'annotation' ARETE entry with subsetting enabled,
with a core similarity threshold of 99% and an accessory similarity of 95%.

```
nextflow run beiko-lab/ARETE \
  --input_sample_table samplesheet.csv \
  --enable_subsetting \
  --core_similarity 99 \
  --accessory_similarity 95 \
  -profile docker \
  -entry annotation
```

Be sure to not include `--skip_poppunk` in your command, because that will then disable all PopPUNK-related processes, including the subsetting subworkflow.
