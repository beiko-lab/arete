# ARETE and dataset size

Currently ARETE has three distinct profiles that change the pipeline execution in some ways: The default profile (which we can call `small`), the `medium` profile and the `large` profile.

These three profiles were developed based on the size of the input dataset and change some parameter defaults based on tests we have performed on similar-sized datasets.

The sizes are:

- For the default or `small` profile, we expect datasets with 100 samples/assemblies or fewer.
  It runs on the default pipeline parameters, with no changes.

- For the `medium` profile, we expect datasets with >100 and <1000 samples.
  It increases the default resource requirements for most processes and also uses [PPanGGoLiN](https://github.com/labgem/PPanGGOLiN) for pangenome construction, instead of [Panaroo](https://github.com/gtonkinhill/panaroo/).

- For the `large` profile, we expect datasets with >1000 samples.
  It also increases default resource requirements for some processes and uses PPanGGoLin.
  Additionally, it enables [PopPUNK subsampling](subsampling.md), with default parameters.
