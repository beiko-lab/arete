process CREATE_REPORT {
    label "process_medium"

    conda (params.enable_conda ? "conda-forge::pandas=1.4.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:ebca4356a18677aaa2c50f396a408343200e514b-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-1e9d4f78feac0eb2c8d8246367973b3f6358defc:ebca4356a18677aaa2c50f396a408343200e514b-0"
    }

    input:
    path annotation
    path diamond_results
    path rgi_output
    path vfdb_fasta
    path phispy_output
    path mobsuite_output

    output:
    path("annotation_report.tsv.gz"), emit: report
    path("feature_profile.tsv.gz"), emit: profile

    script:
    """
    create_report.py \\
        --annotation_out $annotation \\
        --diamond_outs $diamond_results \\
        --rgi_out $rgi_output \\
        --vfdb_fasta $vfdb_fasta \\
        --phispy_out $phispy_output \\
        --mobsuite_out $mobsuite_output
    """

    stub:
    """
    touch annotation_report.tsv.gz
    touch feature_profile.tsv.gz
    """
}
