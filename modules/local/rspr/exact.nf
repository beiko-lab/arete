process RSPR_EXACT {
    label 'process_medium'

    conda "bioconda::ete3=3.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/jvfe/rspr:v1.3.0':
        'docker.io/jvfe/rspr:v1.3.0' }"

    input:
    path subset_df
    path rooted_gene_trees

    output:
    path "exact_output_*csv", emit: csv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    rspr_exact.py \\
        --dataframe $subset_df \\
        $args
    """
    stub:
    """
    touch exact_output_group_1.csv
    """
}
