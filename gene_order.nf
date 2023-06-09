include { EXTRACTION } from './modules/local/gene_order/extraction'

workflow {

    // take:
    // assemblies
    // genbanks
    // rgi_outs

    // main:

    ch_versions = Channel.empty()

    // Optional extraction params
    num_neighbors = params.num_neighbors
    percent_cutoff = params.gene_order_percent_cutoff

    EXTRACTION (
        file(params.input_file_path),
        file(params.rgi_out),
        file(params.gbk_path),
        file(params.gene_order_html_template),
        num_neighbors,
        percent_cutoff
    )


    // emit:

    // versions = ch_versions                     // channel: [ versions.yml ]
}

