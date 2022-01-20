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
//WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { ARETE } from './workflows/arete'

//
// WORKFLOW: Run main nf-core/arete analysis pipeline
//
/*
workflow NFCORE_ARETE {
    ARETE ()
}
*/
workflow assembly{
    include { ASSEMBLY } from './workflows/arete'
    ASSEMBLY ()

}

workflow annotation {
    include { ANNOTATION } from './workflows/arete'
    ANNOTATION ()
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
