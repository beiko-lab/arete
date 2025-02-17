process POPPUNK_QCDB {
    label 'process_high'

    conda (params.enable_conda ? "bioconda::poppunk=2.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.6.0--py39h9b916c0_0':
        'quay.io/biocontainers/poppunk:2.6.0--py39h9b916c0_0' }"

    input:

    path poppunk_db
    val  type_isolate

    output:

    path "poppunk_db", emit: poppunk_db
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
