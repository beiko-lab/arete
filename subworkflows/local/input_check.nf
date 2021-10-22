/*
 * Check input samplesheet and get read channels
 */

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )
include { ASSEMBLYSHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    
    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .splitCsv ( header:true, sep:',' )
        .map { get_sample_info(it) }
        .set { reads }

    emit:
    reads // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def get_sample_info(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return array    
}

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

    emit:
    genomes // channel: [ val(meta), [ reads ] ]
}
// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def get_sample_info_assemblies(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample

    def array = []
    if (!file(row.path).exists()) {
        print("***")
        print(row)
        print("***")
        exit 1, "ERROR: Please check input samplesheet -> Sequence file does not exist!\n${row.fastq_1}"
    }
    array = [ meta, [file(row.path)]]
    
    return array    
}