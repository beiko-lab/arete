process RSPR_EXACT {
    label 'process_single'
    label 'process_long'

    conda "bioconda::ete3=3.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/jvfe/rspr:v1.3.7':
        'docker.io/jvfe/rspr:v1.3.7' }"

    input:
    path subset_df
    path rooted_reference
    path rooted_gene_trees
    val min_branch_length
    val max_support_threshold
    val max_approx_rspr

    output:
    path "exact_output_*tsv", emit: tsv
    path "cluster_file_*txt", emit: txt

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    rspr_exact.py \\
        --dataframe $subset_df \\
        --min_branch_length $min_branch_length \\
        --max_support_threshold $max_support_threshold \\
        --max_approx_rspr $max_approx_rspr \\
        $args
    """
    stub:
    """
    touch exact_output_group_1.tsv
    """
}
