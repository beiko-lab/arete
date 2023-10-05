process RSPR_APPROX {
    label 'process_low'

    conda "bioconda::ete3=3.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/jvfe/rspr:v1.3.0':
        'docker.io/jvfe/rspr:v1.3.0' }"

    input:
    path core_tree
    path gene_tree_list
    val min_rspr_distance

    output:
    path "approx", emit: res_dir
    path "approx/output.csv", emit: csv
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
        --acc \$(cat $gene_tree_list) \\
        -o approx \\
        --min_rspr_distance $min_rspr_distance \\
        $args
    """
    stub:
    """
    mkdir approx
    touch approx/output.csv
    mkdir approx/rooted_gene_trees
    mkdir approx/rooted_reference_tree
    """
}
