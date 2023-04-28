process ADD_GENOME_COLUMN {
    tag "$meta.id"
    label "process_low"

    conda (params.enable_conda ? "conda-forge::pandas=1.4.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pandas:1.4.3"
    } else {
        container "quay.io/biocontainers/pandas:1.4.3"
    }

    input:
    tuple val(meta), path(tsv)
    val dbname
    val skip_n_rows

    output:
    tuple val(meta), path("${prefix}.txt"), emit: txt

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_${dbname}"

    """
    add_column.py $tsv $meta.id $skip_n_rows ${prefix}.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_${dbname}"
    """
    touch ${prefix}.txt
    """
}
