# ARETE Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2 dev - January 19,2021
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
