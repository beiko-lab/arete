process RSPR_HEATMAP {
    label 'process_single'

    conda "bioconda::ete3=3.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/jvfe/rspr:v1.3.4':
        'docker.io/jvfe/rspr:v1.3.4' }"

    input:
    path df

    output:
    path "exact_output.png", emit: png

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    rspr_heatmap.py \\
        --dataframe $df \\
        -o exact_output.png \\
        $args
    """
    stub:
    """
    touch exact_output.png
    """
}
