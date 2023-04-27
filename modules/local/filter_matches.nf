process FILTER_MATCHES {
    tag "$meta.id"
    label "process_low"

    conda (params.enable_conda ? "conda-forge::pandas=1.4.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pandas:1.4.3"
    } else {
        container "quay.io/biocontainers/pandas:1.4.3"
    }

    input:
    tuple val(meta), path(aln)
    val dbname
    val header
    val pident
    val qcover

    output:
    tuple val(meta), path("*.txt"), emit: txt

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${dbname}_filtered"

    """
    filter_alignment.py ${aln} '${meta.id}' '${header}' \\
        ${pident} ${qcover} ${prefix}.txt
    """

    stub:
    """
    touch ${prefix}.txt
    """
}
