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
    def sheet = samplesheet ? "--samplesheet_file $samplesheet": ''
    def columns = samplesheet_columns ? "--samplesheet_columns $samplesheet_columns": ''
    """
    Feature_Dispersion.py \\
        --output_base FeatureDispersion \\
        --tree_file $core_tree \\
        --feature_file $feature_profile \\
        $sheet \\
        $columns
    """
    stub:
    """
    touch FeatureDispersion.tsv
    touch FeatureDispersion.png
    """
}
