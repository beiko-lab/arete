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

include { ORTHOFINDER } from '../../modules/local/orthofinder' addParams (options: [:])


// Usage pattern from nf-core/rnaseq: Empty dummy file for optional inputs
//ch_dummy_input = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/


workflow ORTHOFINDER_PANGENOME{
    take: proteins


}