process CREATE_REPORT {
    label "process_medium"

    conda (params.enable_conda ? "conda-forge::pandas=1.4.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pandas:1.4.3"
    } else {
        container "quay.io/biocontainers/pandas:1.4.3"
    }

    input:
    path annotation
    path diamond_results
    path rgi_output
    path mobsuite_output

    output:
    path("annotation_report.tsv"), emit: tsv

    script:
    """
    create_report.py $annotation $diamond_results \\
        $rgi_output $mobsuite_output
    """

    stub:
    """
    touch annotation_report.tsv
    """
}
