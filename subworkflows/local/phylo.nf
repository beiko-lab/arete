//
// MODULE: Installed directly from nf-core/modules
//
include { IQTREE } from '../../modules/nf-core/iqtree/main'
include { PANAROO_RUN } from '../../modules/nf-core/panaroo/run/main'
include { SNPSITES } from '../../modules/nf-core/snpsites/main'
include { FASTTREE } from '../../modules/nf-core/fasttree/main'

//
// MODULE: Local to the pipeline
//
include { GML2GV } from '../../modules/local/graphviz/gml2gv/main'
include { GET_SOFTWARE_VERSIONS } from '../../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )


workflow PHYLOGENOMICS{
    take:
        gffs
        use_full_alignment
        use_fasttree
    main:
        ch_software_versions = Channel.empty()

        /*
        * Core gene identification and alignment
        */
        // By default, run panaroo
        PANAROO_RUN(gffs.collect{ meta, gff -> gff}.map( gff -> [[id: 'panaroo'], gff]))
        PANAROO_RUN.out.aln.collect{meta, aln -> aln }.set{ ch_core_gene_alignment }
        PANAROO_RUN.out.graph_gml.map{ id, path -> path }.set { panaroo_graph }
        ch_software_versions = ch_software_versions.mix(PANAROO_RUN.out.versions.ifEmpty(null))

        GML2GV(panaroo_graph)

        /*
        * Maximum likelihood core gene tree. Uses SNPSites by default
        */
        if (use_fasttree){
            FASTTREE(ch_core_gene_alignment)
            ch_software_versions = ch_software_versions.mix(FASTTREE.out.versions.ifEmpty(null))
        }

        else{
            if (!use_full_alignment){
                SNPSITES(ch_core_gene_alignment)
                ch_software_versions = ch_software_versions.mix(SNPSITES.out.versions.ifEmpty(null))
                IQTREE(SNPSITES.out.fasta, SNPSITES.out.constant_sites_string)
            }
            else{
                IQTREE(ch_core_gene_alignment, false)
            }
            ch_software_versions = ch_software_versions.mix(IQTREE.out.versions.ifEmpty(null))
        }

    emit:
        phylo_software = ch_software_versions
        core_alignment = ch_core_gene_alignment
}
