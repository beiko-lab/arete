include { PHYLOSHEET_CHECK } from '../../modules/local/samplesheet_check'

// Input to the phylo work flow is different
// Instead of reads, pass in GFF files
workflow PHYLO_INPUT_CHECK {
    take:
    samplesheet

    main:
    PHYLOSHEET_CHECK ( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .map { get_sample_info_phylo(it) }
        .set { genomes }

    //Check that no dots "." are in sample ID
    genomes
        .map { meta, reads -> meta.id }
        .subscribe { if ( "$it".contains(".") ) exit 1, "Please review data input, sampleIDs may not contain dots, but \"$it\" does." }

    emit:
    genomes // channel: [ val(meta), [ reads ] ]
}
// Function to get list of [ meta, [ path ] ]
def get_sample_info_phylo(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end = true

    def array = []
    if (!file(row.path).exists()) {
        print("***")
        print(row)
        print("***")
        exit 1, "ERROR: Please check input samplesheet -> Sequence file does not exist!\n${row.path}"
    }
    array = [ meta, file(row.path)]

    return array
}
