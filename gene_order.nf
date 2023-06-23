include { EXTRACTION } from './modules/local/gene_order/extraction'
include { CLUSTERING } from './modules/local/gene_order/clustering'

include { DIAMOND_MAKEDB } from './modules/nf-core/diamond/makedb/main'
include { DIAMOND_BLASTP as DIAMOND_GENE_ORDER } from './modules/nf-core/diamond/blastp/main'

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

    assembly_channel
        .combine(dbs)
        .set { all_combs }

    DIAMOND_GENE_ORDER (
        all_combs.map{meta, fasta, db -> tuple(meta, fasta)},
        all_combs.map{meta, fasta, db -> db},
        "txt",
        blast_columns
    )

    DIAMOND_GENE_ORDER.out.txt.collect{id, path -> path}.set { blastFiles }

    CLUSTERING (
        assemblyFiles.collect(),
        EXTRACTION.out.fasta_path,
        blastFiles,
        EXTRACTION.out.json_path,
        EXTRACTION.out.neighborhood_indices,
        num_neighbors,
    	inflation,
    	epsilon,
        minpts,
        file(params.gene_order_html_template),
        params.plot_clustering
    )
    // emit:

    // versions = ch_versions                     // channel: [ versions.yml ]
}

