process POPPUNK_MAKE_SAMPLESHEET {
    label "process_low"

    input:
    path(samplesheets)

    output:
    path("poppunk_samplesheet.tsv"), emit: full_samplesheet

    script:
    """
    cat ${samplesheets} > poppunk_samplesheet.tsv
    """

    stub:
    """
    touch poppunk_samplesheet.tsv
    """
}


process POPPUNK_MAKEDB {
    label 'process_high'

    conda (params.enable_conda ? "bioconda::poppunk=2.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.4.0--py39h8884e85_2':
        'quay.io/biocontainers/poppunk:2.4.0--py38hf6d6cf9_1' }"

    input:

    path(filesheet)

    output:

    path("popdb"), emit: poppunk_db
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    poppunk \\
        --create-db \\
        --output popdb \\
        --r-files $filesheet  \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        poppunk: \$(echo \$(poppunk --version 2>&1) | sed 's/^.*poppunk //;')
    END_VERSIONS
    """
}

process POPPUNK_FITMODEL {
    label 'process_high'

    conda (params.enable_conda ? "bioconda::poppunk=2.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.4.0--py39h8884e85_2':
        'quay.io/biocontainers/poppunk:2.4.0--py38hf6d6cf9_1' }"

    input:

    path(poppunk_db)

    output:

    path("popdb_dbscan"), emit: popdb_dbscan
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    poppunk \\
        --fit-model bgmm \\
        --ref-db $poppunk_db \\
        --output popdb_dbscan \\
        --max-a-dist 0.9 \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        poppunk: \$(echo \$(poppunk --version 2>&1) | sed 's/^.*poppunk //;')
    END_VERSIONS
    """
}
