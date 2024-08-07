include { EXTRACTION } from '../../modules/local/gene_order/extraction'
include { CLUSTERING } from '../../modules/local/gene_order/clustering'

include { DIAMOND_MAKEDB } from '../../modules/nf-core/diamond/makedb/main'
include { DIAMOND_BLASTP as DIAMOND_GENE_ORDER } from '../../modules/nf-core/diamond/blastp/main'

workflow GENE_ORDER {

    take:
    assemblies
    genbanks
    rgi_outs

    main:

    rgiFiles = rgi_outs.map { meta, path -> path }.collect()
    gbkFiles = genbanks.map { meta, path -> path }.collect()

    // Optional extraction params
    num_neighbors = params.num_neighbors
    percent_cutoff = params.gene_order_percent_cutoff
    label_cols = params.gene_order_label_cols ? params.gene_order_label_cols : 'None'

    // Optional clustering params
    inflation = params.inflation
    epsilon = params.epsilon
    minpts = params.minpts

    EXTRACTION (
        file(params.input_file_path),
        rgiFiles,
        gbkFiles,
        num_neighbors,
        percent_cutoff,
        label_cols
    )

    // Keep only columns required for clustering module
    def blast_columns = "qseqid sseqid pident bitscore"

    DIAMOND_MAKEDB (
        assemblies.map { meta, path -> path }
    )

    DIAMOND_MAKEDB.out.db.set { dbs }

    // Get all possible combinations between assemblies
    assemblies
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
        assemblies.map { meta, path -> path }.collect(),
        EXTRACTION.out.fasta_path,
        blastFiles,
        EXTRACTION.out.json_path,
        EXTRACTION.out.neighborhood_indices,
        num_neighbors,
    	inflation,
    	epsilon,
        minpts,
        params.plot_clustering
    )

    emit:
    clustering_directory = CLUSTERING.out.cluster_path
    clustering_summary = CLUSTERING.out.summary
}

