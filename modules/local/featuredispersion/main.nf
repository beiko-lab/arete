process FEATURE_DISPERSION {
    label 'process_single'

    conda "bioconda::ete3=3.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/jvfe/rspr:v1.3.7':
        'docker.io/jvfe/rspr:v1.3.7' }"

    input:
    path core_tree
    path feature_profile
    path samplesheet
    val samplesheet_columns

    output:
    path "FeatureDispersion.tsv", emit: tsv
    path "FeatureDispersion.png", emit: png

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Feature_Dispersion.py \\
        --tree_file $core_tree \\
        --feature_file $feature_profile \\
        --samplesheet_file $samplesheet \\
        --samplesheet_columns $samplesheet_columns
    """
    stub:
    """
    touch FeatureDispersion.tsv
    touch FeatureDispersion.png
    """
}
