process VERTICALL_REPAIR {
    tag "$cluster"
    label 'process_high'

    conda "bioconda::verticall=0.4.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/verticall:0.4.1--pyhdfd78af_0':
        'quay.io/biocontainers/verticall:0.4.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${prefix}_repaired.fna"), emit: repaired
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    verticall \\
        repair \\
        $args \\
        -i $assembly \\
        -o ${prefix}_repaired.fna

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verticall: \$(echo \$(verticall --version 2>&1) | sed "s/^Verticall v//g")
    END_VERSIONS
    """
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_repaired.fna

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verticall: \$(echo \$(verticall --version 2>&1) | sed "s/^Verticall v//g")
    END_VERSIONS
    """
}
