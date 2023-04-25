process CONCAT_ALIGNMENT {
    label "process_low"

    conda (params.enable_conda ? "conda-forge::coreutils=8.31" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/coreutils:8.31--h14c3975_0"
    } else {
        container "quay.io/biocontainers/coreutils:8.31--h14c3975_0"
    }

    input:
    path(aln)
    val dbname

    output:
    path("*.txt"), emit: txt

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${dbname}"

    """
    sed -s 1d ${aln} > no_header.txt
    sed -sn 1p ${aln} | uniq > header.txt
    cat header.txt no_header.txt > ${prefix}.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${dbname}"
    """
    touch ${prefix}.txt
    """
}
