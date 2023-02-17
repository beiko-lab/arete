#!/usr/bin/env nextflow
/*
========================================================================================
    fmaguire/arete
========================================================================================
    Github : https://github.com/fmaguire/arete
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

//params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/
// TODO new function to validate our params
WorkflowArete.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { ARETE } from './workflows/arete'
include { ASSEMBLY } from './workflows/arete'
include { ANNOTATION } from './workflows/arete'
include { QUALITYCHECK } from './workflows/arete'


//
// WORKFLOW: Run main nf-core/arete analysis pipeline
//
/*
workflow NFCORE_ARETE {
    ARETE ()
}
*/
workflow assembly{
    ASSEMBLY ()

}

workflow annotation {
    ANNOTATION ()
}

workflow assembly_qc {
    QUALITYCHECK()
}
/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    ARETE ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
