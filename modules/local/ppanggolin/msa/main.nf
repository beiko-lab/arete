process PPANGGOLIN_MSA {
    tag "$meta.id"
    label 'process_high'
    label 'process_high_memory'

    conda "bioconda::ppanggolin=2.0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin:2.0.5--py39hf95cd2a_0':
        'quay.io/biocontainers/ppanggolin:2.0.5--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(pangenome)

    output:
    path "ppanggolin_msa/msa_all_protein/*.aln", emit: alignments
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp $pangenome copied_pangenome.h5

    ppanggolin \\
        msa \\
        $args \\
        --cpu $task.cpus \\
        --pangenome copied_pangenome.h5 \\
        --output ppanggolin_msa \\
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
