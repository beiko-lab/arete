process ADD_GENOME_COLUMN {
    label "process_low"

    conda (params.enable_conda ? "conda-forge::pandas=1.4.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "docker://docker.io/biocontainers/seaborn:0.12.2_cv1"
    } else {
        container "docker.io/biocontainers/seaborn:0.12.2_cv1"
    }

    input:
    path tsv
    val dbname
    val skip_n_rows
    val header_line

    output:
    path("${prefix}.txt"), emit: txt

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${dbname}"

    """
    for table in \$(ls $tsv); do

        sampleid=\$(basename \$table .tsv | sed 's/_${dbname}.*//Ig')

        add_column.py \$table "\${sampleid}" $skip_n_rows "\${sampleid}_added.txt"
    done

    sed -s '1,${header_line}d' *_added.txt > no_header.txt
    sed -sn ${header_line}p *_added.txt | uniq > header.txt
    cat header.txt no_header.txt > ${prefix}.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${dbname}"
    """
    touch ${prefix}.txt
    """
}
