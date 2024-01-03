process RSPR_APPROX {
    label 'process_low'

    conda "bioconda::ete3=3.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/jvfe/rspr:v1.3.4':
        'docker.io/jvfe/rspr:v1.3.4' }"

    input:
    path core_tree
    path gene_tree_list
    path annotation
    val min_rspr_distance
    val min_branch_length
    val max_support_threshold

    output:
    path "approx", emit: res_dir
    path "approx/output.tsv", emit: tsv
    path "approx/rooted_gene_trees", emit: rooted_gene_trees
    path "approx/rooted_reference_tree", emit: rooted_reference_tree
    path "*.csv", emit: csvs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    rspr_approx.py \\
        --core $core_tree \\
        --acc $gene_tree_list \\
        --annotation $annotation \\
        -o approx \\
        --min_rspr_distance $min_rspr_distance \\
        --min_branch_length $min_branch_length \\
        --max_support_threshold $max_support_threshold \\
        $args
    """
    stub:
    """
    mkdir approx
    touch approx/output.tsv
    mkdir approx/rooted_gene_trees
    mkdir approx/rooted_reference_tree
    """
}
