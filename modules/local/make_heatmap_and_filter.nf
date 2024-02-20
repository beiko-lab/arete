process MAKE_HEATMAP_AND_FILTER {
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::seaborn=0.12.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/biocontainers/seaborn:0.12.2_cv1':
        'docker.io/biocontainers/seaborn:0.12.2_cv1' }"

    input:
    path distances
    val core
    val accessory

    output:
    path("thresholds_heatmap.pdf"), emit: plot
    path("genome_mapping.csv"), emit: mapping
    path("removed_genomes.txt"), emit: removed_genomes

    script:
    """
    filter_genomes.py ${distances} \\
        thresholds_heatmap.pdf \\
        --core_threshold $core \\
        --accessory_threshold $accessory
    """

    stub:
    """
    touch thresholds_heatmap.pdf
    touch removed_genomes.txt
    touch genome_mappings.tsv
    """
}
