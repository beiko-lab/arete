include { EXTRACTION } from './modules/local/gene_order/extraction'
include { CLUSTERING } from './modules/local/gene_order/clustering'

include { DIAMOND_MAKEDB } from './modules/nf-core/diamond/makedb/main'
include { DIAMOND_BLASTP } from './modules/nf-core/diamond/blastp/main'

workflow {

    // take:
    // assemblies
    // genbanks
    // rgi_outs

    // main:

    ch_versions = Channel.empty()
    assemblyFiles = Channel.fromPath("${params.assembly_path}/*.{fa,faa,fna}")

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

    // Keep only columns required for clustering module
    def blast_columns = "qseqid sseqid pident bitscore"

    DIAMOND_MAKEDB (
        assemblyFiles
    )

    DIAMOND_MAKEDB.out.db.set { dbs }

    assemblyFiles
        .map {
            item ->
                tuple([id: item.getSimpleName()], item)
        }
        .set { assembly_channel }

    DIAMOND_BLASTP (
        assembly_channel,
        dbs,
        "txt",
        blast_columns
    )

    DIAMOND_BLASTP.out.txt.set { blastFiles }


    // emit:

    // versions = ch_versions                     // channel: [ versions.yml ]
}

