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

include { MAFFT } from '../../modules/nf-core/modules/mafft/main' addParams ( options: [:] )
include { BLASTP } from '../../modules/nf-core/modules/diamond/blastp/main' addParams ( options: [:] )
//
// MODULE: Local to the pipeline
//

include { OF_PREPARE_BLAST } from '../../modules/local/orthofinder' addParams (options: [:])
include { ORTHOFINDER_RUN } from '../../modules/local/orthofinder' addParams (options: [:])


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
process 

workflow ORTHOFINDER_PANGENOME{
    take: proteins

    main:
        ch_software_versions = Channel.empty()

        OF_PREPARE_BLAST(proteins)



}

//function to get list of blast pairs
