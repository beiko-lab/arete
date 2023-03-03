process POPPUNK_MAKE_SAMPLESHEET {
    label "process_low"

    input:
    path(samplesheets)

    output:
    path("poppunk_samplesheet.tsv"), emit: full_samplesheet

    script:
    """
    cat ${samplesheets} > poppunk_samplesheet.tsv
    """

    stub:
    """
    touch poppunk_samplesheet.tsv
    """
}
