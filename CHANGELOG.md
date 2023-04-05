# ARETE Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0 dev - April 5, 2023

### `Added`

- Adds stubs to every process (2a990d48d1d024)
- Adds assembly as a subworkflow (9e6d2679242)
- Adds dbcache (50193b35)
- Updates documentation
  - README, ROADMAP, CITATIONS, docs/usage.md, docs/output.md
- Adds test data (c57297ad922)

- Reorganizes output directories by subworkflow (ff53e97129effb)
- Adds option to skip Kraken2 (192ab2918997bb0d35d2195b8d3da)
- **Adds unit tests to ARETE** (6a10be029268eb2f1d526913a8022a)
- **Adds CI** (49eb2929f287643446c5eb073670435e45a559a9)
- Adds initialization script to the workflow (300f52429ee9cee5e056)
- **Make FastTree the default under the Phylo subwf** (b3a3ffa0e64367)
- Adds light parameter to the workflow (b3a3ffa0e64367)
- Concatenate and filter alignment results (8e738ef237bcc80b6d9)
- **Adds PopPUNK** (d4ea29eaaad54905a81ca5a0c31a699c9ff89b2)
- Updates and improves Bakta, making it the default(3cb354c9c6e20346e6a844e9a0c55060fb8dff)
- Adds DIAMOND alignment against ICEberg2 (ed548cb8ac8bf44d68b)
- **Adds IslandPath** (0c2ba57999f4952e1eb1db3b3107a4ccbdc5ad86)
- **Adds VIBRANT** (f0b1f99a2c7ded75e6684e8775e052ca2ec6e140)
- **Adds Integron Finder** (00735d15aa5d6fef67793fd32152ea91cfa65c2a)
- **Adds Genome subsampling using PopPUNK distances** (af41029c9ce6d2308f3)

### `Fixed`

- Makes assembly and assemblyqc entries available (4fd48439be9)
- Change databases to value channels so they can be recycled (1285c0582b5f78)
- Adds missing params to config and schema (a2c2ea7d216b15f3652)
- Integrate phylo subworkflow to the main pipeline (069400dde1d098fe)

### `Deprecated`

- Removed roary (786c69c223261846)
- Removed orthofinder (8ce706a7a040d0)

## v.1.2.1 dev - March 18, 2022

### `Added`

- Add CheckM nf-core module
- Add `assembly_qc` workflow; runs QC software on existing assemblies

### `Fixed`

- Fixed MultiQC not finding QUAST report (again)

## v1.2.0 dev - January 19,2022

Major refactoring to align with Nextflow DSL2 syntax.

### `Added`

- Add SNPsites module
- Add Diamond nf-core module
- Add IQTree nf-core module

### `Fixed`

- Fixed broken `assembly` workflow entry
- Update most nf-core modules
- Local software modules are now found in scripts rather than directories, as per nf-core v.2.0+ convention
- Update directory structures, manifest files, Groovy files to nf-core v.2.0+ versions
- Update default `max_retries` in config file to 2 (from 1)

### `Dependencies`

- Bump nf-core to version 2.0+. (Note: nf-core is not required to run the pipeline, but module updates, linting etc require it)

### `Deprecated`

- Local diamond module
- Local IQTree module
- Pathracer module currently not used. May be added back later

## v1.1dev - November 3, 2021

Update of tools used with some debugging.

### `Added`

- Add Roary for core genome alignment
- Add Diamond functions for BLASTing VFDB, BacMet
- Add 'annotation only' workflow
- Add 'assembly only' workflow

### `Fixed`

- Fixed MobSuite container
- Fixed IQTree to work with Roary instead of Snippy
- Fixed MultiQC to properly display reports from QUAST, Prokka
- Update singularity containers where relevant

### `Dependencies`

- No changes in dependencies

### `Deprecated`

- Snippy module
- Abricate module
- CRISPR/CasFinder module

## v1.0dev - [date]

Initial release of nf-core/arete, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
