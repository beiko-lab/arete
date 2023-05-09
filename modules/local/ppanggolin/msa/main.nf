process PPANGGOLIN_MSA {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::ppanggolin=1.2.105"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin:1.2.105--py37h8902056_0':
        'biocontainers/ppanggolin:1.2.105--py37h8902056_0' }"

    input:
    tuple val(meta), path(pangenome)

    output:
    tuple val(meta), path("${prefix}_msa")    , emit: results
    path "${prefix}_msa/msa_all_protein/*.aln", emit: alignments
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ppanggolin \\
        msa \\
        $args \\
        --cpu $task.cpus \\
        --pangenome $pangenome \\
        --output "${prefix}"_msa \\
        --partition all

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ppanggolin: \$(echo \$(ppanggolin --version 2>&1) | sed 's/^.*ppanggolin //' ))
    END_VERSIONS
    """
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_msa/msa_all_protein/
    touch ${prefix}_msa/msa_all_protein/sample.aln

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ppanggolin: \$(echo \$(ppanggolin --version 2>&1) | sed 's/^.*ppanggolin //' ))
    END_VERSIONS
    """
}
