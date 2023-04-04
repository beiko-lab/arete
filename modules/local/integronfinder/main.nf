process INTEGRON_FINDER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::integron_finder=2.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/integron_finder:2.0.2--pyhdfd78af_0':
        'quay.io/biocontainers/integron_finder:2.0.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("Results_Integron_Finder_*"), emit: results
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    integron_finder \\
        $args \\
        --cpu $task.cpus \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        integron_finder: \$(echo \$(integron_finder --version 2>&1) | sed -n 's/integron_finder version \\([0-9.]\\+\\).*/\\1/p')
    END_VERSIONS
    """
    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir Results_Integron_Finder_${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        integron_finder: \$(echo \$(integron_finder --version 2>&1) | sed -n 's/integron_finder version \\([0-9.]\\+\\).*/\\1/p')
    END_VERSIONS
    """
}
