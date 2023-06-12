process CONCAT_ALIGNMENT {
    label "process_low"

    conda (params.enable_conda ? "conda-forge::pandas=1.4.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:ebca4356a18677aaa2c50f396a408343200e514b-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:ebca4356a18677aaa2c50f396a408343200e514b-0"
    }

    input:
    path alignments
    path exact_core
    path soft_core

    output:
    path("core_gene_alignment.aln"), emit: core_aln

    script:
    """
    cat $exact_core $soft_core | sort | uniq > full_core.txt
    cat full_core.txt | concatenate_aln.py > core_gene_alignment.aln
    """

    stub:
    """
    touch core_gene_alignment.aln
    """
}
