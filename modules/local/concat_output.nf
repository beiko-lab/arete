process CONCAT_OUTPUT {
    label "process_low"

    conda (params.enable_conda ? "conda-forge::sed=4.7.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/sed:4.7.0"
    } else {
        container "quay.io/biocontainers/sed:4.7.0"
    }

    input:
    path(ann_out)
    val dbname
    val header_line

    output:
    path("*.txt"), emit: txt

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${dbname}"

    """
    sed -s '1,${header_line}d' ${ann_out} > no_header.txt
    sed -sn ${header_line}p ${ann_out} | uniq > header.txt
    cat header.txt no_header.txt > ${prefix}.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${dbname}"
    """
    touch ${prefix}.txt
    """
}
