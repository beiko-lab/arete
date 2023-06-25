process FILTER_AND_CONCAT_MATCHES {
    label "process_low"

    conda (params.enable_conda ? "conda-forge::pandas=1.4.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pandas:1.4.3"
    } else {
        container "quay.io/biocontainers/pandas:1.4.3"
    }

    input:
    path aln
    val dbname
    val header
    val pident
    val qcover
    val header_line

    output:
    path("${prefix}.txt"), emit: txt

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${dbname}"

    """
    for aln in \$(ls $aln); do

        sampleid=\$(basename \$aln .txt | sed 's/_[[:upper:]]*//g')

        filter_alignment.py \$aln '\$sampleid' '${header}' \\
            ${pident} ${qcover} "\${sampleid}_filtered.txt"
    done

    sed -s '1,${header_line}d' *_filtered.txt > no_header.txt
    sed -sn ${header_line}p *_filtered.txt | uniq > header.txt
    cat header.txt no_header.txt > ${prefix}.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${dbname}"
    """
    touch ${prefix}.txt
    """
}
