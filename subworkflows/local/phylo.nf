//
// MODULE: Installed directly from nf-core/modules
//
include { IQTREE } from '../../modules/nf-core/iqtree/main'
include { PANAROO_RUN } from '../../modules/nf-core/panaroo/run/main'
include { SNPSITES } from '../../modules/nf-core/snpsites/main'
include { FASTTREE as CORE_FASTTREE } from '../../modules/nf-core/fasttree/main'

//
// MODULE: Local to the pipeline
//
include { CHUNKED_FASTTREE as GENE_FASTTREE } from '../../modules/local/chunked_fasttree'
include { PPANGGOLIN_WORKFLOW } from '../../modules/local/ppanggolin/workflow/main'
include { PPANGGOLIN_MSA } from '../../modules/local/ppanggolin/msa/main'
include { GML2GV } from '../../modules/local/graphviz/gml2gv/main'
include { GET_SOFTWARE_VERSIONS } from '../../modules/local/get_software_versions'
include { CONCAT_ALIGNMENT } from '../../modules/local/concat_alignment'


workflow PHYLOGENOMICS{
    take:
        gffs
        use_full_alignment
        use_fasttree
    main:
        ch_software_versions = Channel.empty()

        gffs
            .map { meta, path -> [meta.id, path.getName()] }
            .set { gff_paths }

        // Create samplesheet
        if (params.use_ppanggolin) {
            gff_paths
                .collectFile(newLine: true) { item ->
                    [ "${item[0]}.txt", item[0] + '\t' + item[1] ]
                }
                .collectFile(name: 'ppanggolin_samplesheet.tsv')
                .map{ path -> [[id: 'ppanggolin'], path] }
                .set { gff_samplesheet }

            PPANGGOLIN_WORKFLOW(gff_samplesheet, gffs.collect{ it[1] })

            PPANGGOLIN_WORKFLOW.out.pangenome.set { pangenome }

            PPANGGOLIN_MSA(pangenome)

            PPANGGOLIN_MSA.out.alignments.set{ ch_all_alignments }

            ch_all_alignments
                .collect()
                .flatten()
                .map{it -> it.toString() }
                .collectFile(name: 'aln_paths.txt', newLine: true)
                .set{ alignment_sheet }

            CONCAT_ALIGNMENT (
                alignment_sheet,
                PPANGGOLIN_WORKFLOW.out.exact_core,
                PPANGGOLIN_WORKFLOW.out.soft_core
            )

            CONCAT_ALIGNMENT.out.core_aln.set { ch_core_alignment }

            ch_software_versions = ch_software_versions.mix(PPANGGOLIN_MSA.out.versions.ifEmpty(null))
        } else {
            gff_paths
                .collectFile(newLine: true) { item ->
                    [ "${item[0]}.txt", item[1] ]
                }
                .collectFile(name: 'panaroo_samplesheet.txt')
                .map{ path -> [[id: 'panaroo'], path] }
                .set { gff_samplesheet }

            /*
            * Core gene identification and alignment
            */
            // By default, run panaroo
            PANAROO_RUN(gff_samplesheet, gffs.collect{ it[1] })
            PANAROO_RUN.out.aln.collect{ meta, aln -> aln }.set{ ch_core_alignment }
            PANAROO_RUN.out.accessory_aln.set{ ch_all_alignments }

            PANAROO_RUN.out.graph_gml.map{ id, path -> path }.set { panaroo_graph }
            ch_software_versions = ch_software_versions.mix(PANAROO_RUN.out.versions.ifEmpty(null))

            GML2GV(panaroo_graph)
        }

        ch_all_alignments
            .flatten()
            .buffer( size: 300, remainder: true )
            .set { chunked_alignments }
        /*
        * Maximum likelihood core gene tree.
        */
        gene_trees = Channel.empty()
        if (use_fasttree){
            def is_nt = params.use_ppanggolin ? false : true
            GENE_FASTTREE(chunked_alignments, is_nt)
            GENE_FASTTREE.out.phylogeny.set { gene_trees }
            ch_software_versions = ch_software_versions.mix(GENE_FASTTREE.out.versions.ifEmpty(null))

            CORE_FASTTREE(ch_core_alignment, is_nt)
            CORE_FASTTREE.out.phylogeny.set { core_tree }
        }

        if (!use_fasttree && !params.use_ppanggolin) {
            if (!use_full_alignment){
                SNPSITES(ch_core_alignment)
                ch_software_versions = ch_software_versions.mix(SNPSITES.out.versions.ifEmpty(null))
                IQTREE(SNPSITES.out.fasta, SNPSITES.out.constant_sites_string)
            }
            else{
                IQTREE(ch_core_alignment, false)
            }
            core_tree = IQTREE.out.phylogeny
            ch_software_versions = ch_software_versions.mix(IQTREE.out.versions.ifEmpty(null))
        }

    emit:
        phylo_software = ch_software_versions
        all_alignments = ch_all_alignments
        core_tree = core_tree
        gene_trees = gene_trees
}
