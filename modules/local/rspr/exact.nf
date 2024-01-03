process RSPR_EXACT {
    label 'process_low'
    label 'process_long'

    conda "bioconda::ete3=3.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/jvfe/rspr:v1.3.2':
        'docker.io/jvfe/rspr:v1.3.2' }"

    input:
    path subset_df
    path rooted_gene_trees
    val min_branch_length
    val max_support_threshold
    val max_approx_rspr

    output:
    path "exact_output_*tsv", emit: tsv

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
