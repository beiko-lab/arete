process POPPUNK_FITMODEL {
    label 'process_high'

    conda (params.enable_conda ? "bioconda::poppunk=2.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.6.0--py39h9b916c0_0':
        'quay.io/biocontainers/poppunk:2.6.0--py39h9b916c0_0' }"

    input:

    path poppunk_db
    val model

    output:

    path "poppunk_${model}" , emit: poppunk_results
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    poppunk \\
        --fit-model $model \\
        --ref-db $poppunk_db \\
        --output poppunk_${model} \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        poppunk: \$(echo \$(poppunk --version 2>&1) | sed 's/^.*poppunk //;')
    END_VERSIONS
    """
}
