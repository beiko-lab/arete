workflow GENEORDER_INPUT_CHECK {
    take:
    samplesheet

    main:
    samplesheet
        .splitCsv(header: true)
        .map { it -> get_sample_info_geneorder(it) }
        .set { geneorder_input }

    geneorder_input
        .map { meta, assemblies, gbks, rgis -> meta.id }
        .subscribe { if ( "$it".contains(".") ) exit 1, "Please review data input, sampleIDs may not contain dots, but \"$it\" does." }


    emit:
    geneorder_input
}

def get_sample_info_geneorder(row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end = true //Bit of a hack; call assemblies "single end" to allow passing to kraken

    def array = []
    if (!file(row.fna_file_path).exists()) {
        print("***")
        print(row.fna_file_path)
        print("***")
        exit 1, "ERROR: Please check input samplesheet -> Assembly file does not exist!\n${row.fna_file_path}"
    }
    if (!file(row.gbk).exists()) {
        print("***")
        print(row.gbk)
        print("***")
        exit 1, "ERROR: Please check input samplesheet -> GenBank file does not exist!\n${row.gbk}"
    }
    if (!file(row.rgi).exists()) {
        print("***")
        print(row.rgi)
        print("***")
        exit 1, "ERROR: Please check input samplesheet -> RGI file does not exist!\n${row.rgi}"
    }
    array = [ meta, file(row.fna_file_path), file(row.gbk), file(row.rgi) ]

    return array
}
