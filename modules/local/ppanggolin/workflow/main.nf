process PPANGGOLIN_WORKFLOW {
    tag "$meta.id"
    label 'process_high'
    label 'process_high_memory'
    label 'process_long'

    conda "bioconda::ppanggolin=2.0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin:2.0.5--py39hf95cd2a_0':
        'quay.io/biocontainers/ppanggolin:2.0.5--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(samplesheet)
    path gffs

    output:
    tuple val(meta), path("$prefix")             , emit: results
    tuple val(meta), path("$prefix/pangenome.h5"), emit: pangenome
    path("$prefix/partitions/exact_core.txt")    , emit: exact_core
    path("$prefix/partitions/soft_core.txt")     , emit: soft_core
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ppanggolin \\
        workflow \\
        $args \\
        --cpu $task.cpus \\
        --anno $samplesheet \\
        --output $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ppanggolin: \$(echo \$(ppanggolin --version 2>&1) | sed 's/^.*ppanggolin //' ))
    END_VERSIONS
    """
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p $prefix/partitions
    touch $prefix/partitions/exact_core.txt
    touch $prefix/partitions/soft_core.txt
    touch $prefix/pangenome.h5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ppanggolin: \$(echo \$(ppanggolin --version 2>&1) | sed 's/^.*ppanggolin //' ))
    END_VERSIONS
    """
}
