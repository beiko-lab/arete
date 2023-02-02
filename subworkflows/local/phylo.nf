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
include { PANAROO_RUN } from '../../modules/nf-core/modules/panaroo/run/main' addParams( options: [args:'-a pan --clean-mode strict'])
include { SNPSITES } from '../../modules/nf-core/modules/snpsites/main' addParams( options: [:] )
include { FASTTREE } from '../../modules/nf-core/modules/fasttree/main' addParams( options: [:] )

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )


workflow PHYLOGENOMICS{
    take:
        gffs
        use_roary
        use_full_alignment
        use_fasttree
    main:
        ch_software_versions = Channel.empty()

        /*
        * Core gene identification and alignment
        */
        if(use_roary){
            ROARY(gffs.collect{ meta, gff -> gff}.map( gff -> [[id: 'roary'], gff]))
            ROARY.out.aln.collect{ meta, aln -> aln }.set{ ch_core_gene_alignment }
            ch_software_versions = ch_software_versions.mix(ROARY.out.versions.ifEmpty(null))
        }

        // By default, run panaroo
        else{
            PANAROO_RUN(gffs.collect{ meta, gff -> gff}.map( gff -> [[id: 'panaroo'], gff]))
            PANAROO_RUN.out.aln.collect{meta, aln -> aln }.set{ ch_core_gene_alignment }
            ch_software_versions = ch_software_versions.mix(PANAROO_RUN.out.versions.ifEmpty(null))
        }

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

/*
workflow MULTI_SPECIES_PHYLO{

}
*/
