#!/usr/bin/env nextflow

params.tree = file("$baseDir/test/rooted_vtec_phylogeny.tre")
params.feat_table = file("$baseDir/test/feature_profile_reallyreduced.tsv")

include { EVOLCCM } from './modules/local/evolccm.nf'

workflow {
    EVOLCCM (
        params.tree,
        params.feat_table
    )
}
