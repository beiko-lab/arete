process ADD_COLUMN {
    tag "$meta.id"
    label "process_low"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    tuple val(meta), path(aln)
    val dbname
    val header

    output:
    tuple val(meta), path("*.txt"), emit: txt

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${dbname}"

    """
    sed "s/\$/\\t${meta.id}/" ${aln} > ${prefix}_addedcol.temp
    echo ${header} > header.temp
    tr ' ' '\\t' < header.temp > header_tabs.temp
    cat header_tabs.temp ${prefix}_addedcol.temp > ${prefix}_addedcol.txt
    """

    stub:
    """
    touch ${prefix}.txt
    """
}
