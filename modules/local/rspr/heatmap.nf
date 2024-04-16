process RSPR_HEATMAP {
    label 'process_single'

    conda "bioconda::ete3=3.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/jvfe/rspr:v1.3.7':
        'docker.io/jvfe/rspr:v1.3.7' }"

    input:
    path df
    path cluster_file
    path rooted_reference
    val min_heatmap_exact_rspr
    val max_heatmap_exact_rspr

    output:
    path "exact_output.png", emit: png
    path "exact_group_output.png", emit: exact_group_output
    path "cluster_tree_output.png", emit: cluster_tree_output
    path "cluster_file_output.nwk", emit: cluster_file_output

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    rspr_heatmap.py \\
        --dataframe $df \\
        --cluster_file $cluster_file \\
        -o exact_output.png \\
        -go exact_group_output.png \\
        -co cluster_tree_output.png \\
        -cfo cluster_file_output.nwk \\
        --min_heatmap_exact_rspr $min_heatmap_exact_rspr \\
        --max_heatmap_exact_rspr $max_heatmap_exact_rspr \\
        $args
    """
    stub:
    """
    touch exact_output.png
    touch exact_group_output.png
    touch cluster_tree_output.png
    touch cluster_file_output.nwk
    """
}
