process MAKE_HEATMAP {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::seaborn=0.12.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/biocontainers/seaborn:0.12.2_cv1':
        'docker.io/biocontainers/seaborn:0.12.2_cv1' }"

    input:
    path distances

    output:
    path("thresholds_heatmap.pdf"), emit: plot

    script:
    """
    make_heatmap.py ${distances} thresholds_heatmap.pdf
    """

    stub:
    """
    touch thresholds_heatmap.pdf
    """
}
