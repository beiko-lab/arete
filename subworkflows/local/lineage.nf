def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { IQTREE } from '../../modules/nf-core/iqtree/main'
include { PANAROO_RUN } from '../../modules/nf-core/panaroo/main'
include { SNPSITES } from '../../modules/nf-core/snpsites/main'
include { FASTTREE } from '../../modules/nf-core/fasttree/main'

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )


workflow SINGLE_SPECIES_PHYLO{
    take:
        gffs
        use_full_alignment
        use_fasttree
    main:
        ch_software_versions = Channel.empty()
        ch_core_gene_alignment = Channel.empty()

}
