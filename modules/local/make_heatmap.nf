process MAKE_HEATMAP {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::vibrant=1.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vibrant:1.2.1--hdfd78af_4':
        'quay.io/biocontainers/vibrant:1.2.1--hdfd78af_4' }"

    input:
    path distances

    output:
    path("thresholds_heatmap.pdf"), emit: plot

    script:
    def args = task.ext.args ?: ''

    """
    make_heatmap.py ${distances} thresholds_heatmap.pdf
    """

    stub:
    """
    touch thresholds_heatmap.pdf
    """
}
