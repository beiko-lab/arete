process BUILD_UPSET {
    label 'process_single'

    conda "conda-forge::upsetplot=0.9.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/jvfe/arete_upset:v0.1.0':
        'docker.io/jvfe/arete_upset:v0.1.0' }"

    input:
    path feature_profile
    path samplesheet
    val samplesheet_columns

    output:
    path "BuildUpset*png", emit: pngs

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Build_Upset.py \\
        --output_file_prefix BuildUpset \\
        --feature_profile_file $feature_profile \\
        --samplesheet_file $samplesheet \\
        --columns $samplesheet_columns
    """
    stub:
    """
    touch BuildUpset_dummy.png
    """
}
