# ARETE Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0 - November 30, 2023

### `Added`

- Automate phylogenomic analysis from ppanggolin outputs (d78dfc9d811bac6f8f0f5d0d71cb31a24dadaf02)
- Add parameter to skip pangenomics/phylogenomics subworkflow (3f3f5b218e94adb6ba6d743e3b8854bebae02eb4)
- Add new ARETE logo (4a3b49275d124ca92186d9e82ea2e700bfbf1f2f)
- **Added `--annotation_tools` parameter** (2cf57a3a06be53d0550a004101d199fc1760f84d)
- Validate sample IDs in annotation input check (09400c0e10919001286c8e829299c03129dd64c2)
- General improvements to documentation (3bdea3083d9dfc3f23a198a23843820ac68be2ca, 21b853bf065da239247ed62a2a670f2d25a56de, a109c3c19576cffa9a1029e58229c3054806a017, a10f7bd4ff73bcadb389d18a82d18c6508792237)
- **Add phylogenomics entry** (1cb656372a752bb9288b3e220cbddc71943d4496)
- **Add recombination subworkflow** (1d10d96612c7c7e7f0a7c703771b1ce7d76b35a9)
- Add QC assembly filtering (ccc2f50c7342bb05598cf8deb64cab332acbcfdb)
- Improve MultiQC reports (b097f9b58722d7ff99a1bb4c057152f4e0875dca)
- **Add rSPR** (b0d0cf1467c7f86aa01e023b359e0fe31e7029ad, 7e494ab18d3bcff2e9da87af103c9b2ebe164b20)
- Add entries for EvolCCM and recombination (8bdad860fd75f6d1b06034201874f4739a20895f, 710ba2524f7850938039d6dd4703b47065d49768)
- **Add Gene Order** (6b7e2685143db2d1153ce8cd58b0baab84b064dc, b22988b7a0b945c63f8261fb27dc4b6ccc0f64b0)

### `Fixed`

- Stop ppaggoolin from altering input file (53cbe771c8af23432cfec09b95be7ac6a27c3d0d)
- Return PopPUNK create db outputs (339620c9257e7d99c8b5f933f298e3b2f3e96b14)
- Fix large dataset error on rSPR entry (9a119b4fa588707e36114d864a0ad4cc38cc239c)

### `Deprecated`

- No longer includes phispy information in the annotation report (82bb3f8069a714cd423bef8bdd9e03b2f0f8f919)

## v0.9 - June 6, 2023

### `Added`

- Make subsampling outputs only be used in Phylogenomics subwf (71240fdf1db73991e9e0ae74cacce3b1007d324b)
- Add gml2gv module (4775118f113e38741852a718304241bfdb57cb60)
- Add entry suggestion (47e4f8589c5cc2370033e6f18b9da6251193cb5b)
- Create concatenated outputs for annotation steps (a15cb80a79311f6ba7a860cc8f716a6f61df1610)
- Add genome similarity heatmap (1dfa297a1c6236f30d83e884bcc5f2426d7cba07)
- Add pident and qcover filtering (e289b5caa624f98e1bcde2ac7d105237ce81c7a2)
- Run RGI on Bakta outputs (de88f6f574b6786eb47bba237894d50baba9886d)
- **Add report with annotation outputs** (3c96481090a0785fbc280a6818be05df2ab33000)
- **Add ppanggolin** (1409df676e9c964c972132f8466840a75b5981d9)
- **Add presence/absence table** (1e43d5032eb026877f18e800ea16bfc8b0821c9d)
- **Add PhiSpy** (080d88cde7cd05e12e40fa1703f4ab77f28e071f)

### `Fixed`

- Make CheckM detect file extensions (4f2067263f84959f96a203e1b18d1d0832d2efa4)
- Make poppunk flag only be checked at relevant points (290df03b07479ec20c077bfc926aadb860aaa96d)
- Change assembly samplesheet message (893c5ace3755e43efafd6ecc075194a7c4eeb614)
- Use panaroo samplesheet (4bcd83371c32ef93be9604c3de116bf4d0867df0)

- ### `Deprecated`

- **Removed VIBRANT** (080d88cde7cd05e12e40fa1703f4ab77f28e071f)

## v2.0.1 dev - April 10, 2023

### `Fixed`

- Makes CheckM auto-detect FASTA file extensions
  - So it doesn't fail with files without `.fna`

### `Added`

- Moves Assembly QC subworkflow into its own file

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
