process GET_RECOMB_INPUT {
    label "process_low"

    conda (params.enable_conda ? "conda-forge::pandas=1.4.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:ebca4356a18677aaa2c50f396a408343200e514b-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:ebca4356a18677aaa2c50f396a408343200e514b-0"
    }

    input:
    path quast
    path poppunk
    path assembly_samplesheet

    output:
    path("recomb_samplesheet.csv"), emit: csv
    path("cluster*txt"), emit: clusters

    script:
    def args = task.ext.args ?: ''

    """
    organize_recomb_data.py \\
        $quast \\
        $poppunk \\
        $assembly_samplesheet \\
        "recomb_samplesheet.csv"
    """

    stub:
    """
    touch recomb_samplesheet.csv
    """
}
