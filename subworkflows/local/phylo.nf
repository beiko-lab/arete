//
// MODULE: Installed directly from nf-core/modules
//
include { IQTREE } from '../../modules/nf-core/iqtree/main'
include { PANAROO_RUN } from '../../modules/nf-core/panaroo/run/main'
include { SNPSITES } from '../../modules/nf-core/snpsites/main'

//
// MODULE: Local to the pipeline
//
include { CHUNKED_FASTTREE as FASTTREE } from '../../modules/local/chunked_fasttree'
include { PPANGGOLIN_WORKFLOW } from '../../modules/local/ppanggolin/workflow/main'
include { PPANGGOLIN_MSA } from '../../modules/local/ppanggolin/msa/main'
include { GML2GV } from '../../modules/local/graphviz/gml2gv/main'
include { GET_SOFTWARE_VERSIONS } from '../../modules/local/get_software_versions'


workflow PHYLOGENOMICS{
    take:
        gffs
        use_full_alignment
        use_fasttree
    main:
        ch_software_versions = Channel.empty()

        gffs
            .map { meta, path -> [meta.id, path.toString()] }
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

            PPANGGOLIN_WORKFLOW(gff_samplesheet)

            PPANGGOLIN_WORKFLOW.out.pangenome.set { pangenome }

            PPANGGOLIN_MSA(pangenome)

            ch_software_versions = ch_software_versions.mix(PPANGGOLIN_MSA.out.versions.ifEmpty(null))
            PPANGGOLIN_MSA.out.alignments.flatten().set{ ppanggolin_msas }
            ppanggolin_msas.buffer( size: 300, remainder: true ).set { ch_gene_alignments }

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
            PANAROO_RUN(gff_samplesheet)
            PANAROO_RUN.out.aln.collect{meta, aln -> aln }.set{ ch_gene_alignments }
            PANAROO_RUN.out.graph_gml.map{ id, path -> path }.set { panaroo_graph }
            ch_software_versions = ch_software_versions.mix(PANAROO_RUN.out.versions.ifEmpty(null))

            GML2GV(panaroo_graph)
        }


        /*
        * Maximum likelihood core gene tree. Uses SNPSites by default
        */
        if (use_fasttree){
            def is_nt = params.use_ppanggolin ? false : true
            FASTTREE(ch_gene_alignments, is_nt)
            ch_software_versions = ch_software_versions.mix(FASTTREE.out.versions.ifEmpty(null))
        }

        if (!use_fasttree && !params.use_ppanggolin) {
            if (!use_full_alignment){
                SNPSITES(ch_gene_alignments)
                ch_software_versions = ch_software_versions.mix(SNPSITES.out.versions.ifEmpty(null))
                IQTREE(SNPSITES.out.fasta, SNPSITES.out.constant_sites_string)
            }
            else{
                IQTREE(ch_gene_alignments, false)
            }
            ch_software_versions = ch_software_versions.mix(IQTREE.out.versions.ifEmpty(null))
        }

    emit:
        phylo_software = ch_software_versions
        core_alignment = ch_gene_alignments
}
