process SKA {
    tag "$cluster"
    label 'process_medium'

    conda "bioconda::gubbins=3.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gubbins:3.3.0--py310pl5321h8472f5a_0':
        'quay.io/biocontainers/gubbins:3.3.0--py310pl5321h8472f5a_0' }"

    input:
    tuple val(cluster), path(assemblies), path(reference)
    path assembly_files

    output:
    tuple val(cluster), path("${cluster}_alignment.aln"), emit: aln
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION = '3.3.0'
    """
    generate_ska_alignment.py \\
        $args \\
        --input $assemblies \\
        --reference $reference \\
        --out ${cluster}_alignment.aln \\
        --threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: $VERSION
    END_VERSIONS
    """
    stub:
    def VERSION = '3.3.0'
    """
    touch ${cluster}_alignment.aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gubbins: $VERSION
    END_VERSIONS
    """
}
