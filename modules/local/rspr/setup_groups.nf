process SETUP_GROUPS {
    label 'process_single'

    conda "bioconda::ete3=3.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/jvfe/rspr:v1.3.0':
        'docker.io/jvfe/rspr:v1.3.0' }"

    input:
    path csv

    output:
    path "rspr_*csv", emit: csvs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    setup_rspr_groups.py \\
        --dataframe $csv
    """
    stub:
    """
    touch rspr_group_1.csv
    """
}
