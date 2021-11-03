# ARETE Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
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
