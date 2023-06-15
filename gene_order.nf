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
    rgiFiles = Channel.fromPath("${params.rgi_out}/*_rgi.txt").collect()
    gbkFiles = Channel.fromPath("${params.gbk_path}/*.{gbk,gbff}").collect()

    // Optional extraction params
    num_neighbors = params.num_neighbors
    percent_cutoff = params.gene_order_percent_cutoff

    // Optional clustering params
    inflation = params.inflation
    epsilon = params.epsilon
    minpts = params.minpts

    EXTRACTION (
        file(params.input_file_path),
        rgiFiles,
        gbkFiles,
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

    DIAMOND_BLASTP.out.txt.collect{id, path -> path}.set { blastFiles }

    // CLUSTERING (
    //     blastFiles,
    //     file(params.assembly_path),
    //     EXTRACTION.out.fasta_path,
    //     num_neighbors,
    // 	inflation,
    // 	epsilon,
    // 	minpts
    // )
    // emit:

    // versions = ch_versions                     // channel: [ versions.yml ]
}

