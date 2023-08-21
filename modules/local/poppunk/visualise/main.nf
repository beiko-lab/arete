process POPPUNK_VISUALISE {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::poppunk=2.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.6.0--py39h9b916c0_0':
        'quay.io/biocontainers/poppunk:2.6.0--py39h9b916c0_0' }"

    input:

    path poppunk_refdb
    path poppunk_querydb

    output:

    path "poppunk_visualizations" , emit: poppunk_visualizations
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    poppunk_visualise \\
        --ref-db $poppunk_refdb \\
        --query-db $poppunk_querydb \\
        --output poppunk_visualizations \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        poppunk: \$(echo \$(poppunk_visualise --version 2>&1) | sed 's/^.*poppunk_visualise //;')
    END_VERSIONS
    """
}
