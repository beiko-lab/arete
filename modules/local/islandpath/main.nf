process ISLANDPATH {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::islandpath=1.0.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/islandpath:1.0.6--hdfd78af_0':
        'quay.io/biocontainers/islandpath:1.0.6--hdfd78af_0' }"

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${prefix}.gff"), emit: gff
    path "Dimob.log"                      , emit: log
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION = '1.0.6'
    """
    islandpath \\
        $genome \\
        ${prefix}.gff \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        islandpath: $VERSION
    END_VERSIONS
    """
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.6'
    """
    touch ${prefix}.gff
    touch Dimob.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        islandpath: $VERSION
    END_VERSIONS
    """
}
