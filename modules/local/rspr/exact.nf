process RSPR_EXACT {
    label 'process_medium'

    conda "bioconda::ete3=3.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/jvfe/rspr:v1.3.0':
        'docker.io/jvfe/rspr:v1.3.0' }"

    input:
    path core_tree
    path gene_tree_list

    output:
    path "approx"

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    rspr_exact.py \\
        --core $core_tree \\
        --acc $gene_tree_list \\
        -o approx \\
        $args
    """
    stub:
    """
    mkdir approx
    """
}
