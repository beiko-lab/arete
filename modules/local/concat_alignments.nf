process CONCAT_ALIGNMENT {
    label "process_low"

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path(aln)
    val dbname

    output:
    path("*.txt"), emit: txt

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${dbname}"

    """
    cat ${aln} > ${prefix}.txt
    """

    stub:
    """
    touch ${prefix}.txt
    """
}
