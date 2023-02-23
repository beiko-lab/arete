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

    conda (params.enable_conda ? "bioconda::poppunk=2.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.6.0--py39h9b916c0_0':
        'quay.io/biocontainers/poppunk:2.6.0--py39h9b916c0_0' }"

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

process POPPUNK_QCDB {
    label 'process_high'

    conda (params.enable_conda ? "bioconda::poppunk=2.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.6.0--py39h9b916c0_0':
        'quay.io/biocontainers/poppunk:2.6.0--py39h9b916c0_0' }"

    input:

    path(poppunk_db)
    val  type_isolate


    output:

    path("popdb"), emit: poppunk_db
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def isolate_arg = type_isolate ? "--type-isolate ${type_isolate}" : ""

    """
    poppunk \\
        --qc-db \\
        --ref-db  $poppunk_db \\
        $isolate_arg \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        poppunk: \$(echo \$(poppunk --version 2>&1) | sed 's/^.*poppunk //;')
    END_VERSIONS
    """
}

process POPPUNK_FITMODEL {
    label 'process_high'

    conda (params.enable_conda ? "bioconda::poppunk=2.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.6.0--py39h9b916c0_0':
        'quay.io/biocontainers/poppunk:2.6.0--py39h9b916c0_0' }"

    input:

    path(poppunk_db)
    val model

    output:

    path("popdb_*"), emit: popdb_dbscan
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "popdb_${model}"

    """
    poppunk \\
        --fit-model $model \\
        --ref-db $poppunk_db \\
        --output ${prefix} \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        poppunk: \$(echo \$(poppunk --version 2>&1) | sed 's/^.*poppunk //;')
    END_VERSIONS
    """
}
