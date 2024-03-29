process SAMPLESHEET_CHECK {
    tag "$samplesheet"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path samplesheet

    output:
    path '*.csv'

    script:
    """
    check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv
    """
}

process ASSEMBLYSHEET_CHECK {
    tag "$assemblysheet"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path assemblysheet

    output:
    path '*.csv'

    script:
    """
    check_assembly_samplesheet.py $assemblysheet assemblysheet.valid.csv
    """
}

process PHYLOSHEET_CHECK {
    tag "$phylosheet"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path phylosheet

    output:
    path '*.csv'

    script:
    """
    check_phylo_samplesheet.py $phylosheet phylosheet.valid.csv
    """
}
