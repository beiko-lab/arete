include { ASSEMBLYSHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )

// Input to the annotation work flow is different
// Instead of reads, pass in genomes or assemblies
workflow ANNOTATION_INPUT_CHECK{
    take:
    samplesheet

    main:
    ASSEMBLYSHEET_CHECK ( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .map { get_sample_info_assemblies(it) }
        .set { genomes }

    //Check that no dots "." are in sample ID
    genomes
        .map { meta, reads -> meta.id }
        .subscribe { if ( "$it".contains(".") ) exit 1, "Please review data input, sampleIDs may not contain dots, but \"$it\" does." }

    emit:
    genomes // channel: [ val(meta), [ reads ] ]
}
// Function to get list of [ meta, [ path ] ]
def get_sample_info_assemblies(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end = true //Bit of a hack; call assemblies "single end" to allow passing to kraken

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
