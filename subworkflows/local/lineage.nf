def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

def modules = params.modules.clone()

//
// MODULE: Installed directly from nf-core/modules
//
include { IQTREE } from '../../modules/nf-core/modules/iqtree/main'  addParams( options: [:] )
include { ROARY } from '../../modules/nf-core/modules/roary/main'  addParams( options: [args:'-e -n'] )
include { PANAROO_RUN } from '../../modules/nf-core/modules/panaroo/main' addParams( options: [args:'-a core'])
include { SNPSITES } from '../../modules/nf-core/modules/snpsites/main' addParams( options: [:] )
include { FASTTREE } from '../../modules/nf-core/modules/fasttree/main' addParams( options: [:] )

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )


workflow SINGLE_SPECIES_PHYLO{
    take:
        gffs
        use_roary
        use_full_alignment
        use_fasttree
    main:
        ch_software_versions = Channel.empty()
        ch_core_gene_alignment = Channel.empty()

        