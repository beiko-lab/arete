
/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input_sample_table, params.multiqc_config, params.reference_genome ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { PHYLO_INPUT_CHECK } from '../subworkflows/local/phylo_input_check'
include { ANNOTATION_INPUT_CHECK } from '../subworkflows/local/annotation_input_check'
include { RSPR_INPUT_CHECK } from '../subworkflows/local/rspr_input_check'
include { GENEORDER_INPUT_CHECK } from '../subworkflows/local/geneorder_input_check'
include { ASSEMBLE_SHORTREADS } from '../subworkflows/local/assembly'
include { ANNOTATE_ASSEMBLIES } from '../subworkflows/local/annotation'
include { CHECK_ASSEMBLIES } from '../subworkflows/local/assemblyqc'
include { PHYLOGENOMICS } from '../subworkflows/local/phylo'
include { RUN_POPPUNK } from '../subworkflows/local/poppunk'
include { RECOMBINATION } from '../subworkflows/local/recombination'
include { SUBSET_GENOMES } from '../subworkflows/local/subsample'
include { EVOLCCM } from '../subworkflows/local/evolccm'
include { RSPR } from '../subworkflows/local/rspr'
include { GENE_ORDER } from '../subworkflows/local/gene_order'
/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC as RAW_FASTQC  } from '../modules/nf-core/fastqc/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { FASTQC as TRIM_FASTQC } from '../modules/nf-core/fastqc/main'
include { FASTP                 } from '../modules/nf-core/fastp/main'
include { UNICYCLER             } from '../modules/nf-core/unicycler/main'
include { PROKKA                } from '../modules/nf-core/prokka/main'
include { GET_CAZYDB;
          GET_VFDB;
          GET_BACMET} from '../modules/local/blast_databases.nf'
include { DIAMOND_MAKEDB as DIAMOND_MAKE_CAZY;
          DIAMOND_MAKEDB as DIAMOND_MAKE_VFDB;
          DIAMOND_MAKEDB as DIAMOND_MAKE_BACMET } from '../modules/nf-core/diamond/makedb/main'
include { DIAMOND_BLASTX as DIAMOND_BLAST_CAZY;
          DIAMOND_BLASTX as DIAMOND_BLAST_VFDB;
          DIAMOND_BLASTX as DIAMOND_BLAST_BACMET } from '../modules/nf-core/diamond/blastx/main'
include { IQTREE } from '../modules/nf-core/iqtree/main'
include { SNPSITES } from '../modules/nf-core/snpsites/main'
//
// MODULE: Local to the pipeline
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { RGI;
          UPDATE_RGI_DB } from '../modules/local/rgi'
include { MOB_RECON } from '../modules/local/mobsuite'
include { KRAKEN2_DB } from '../modules/local/get_minikraken'
include { GET_DB_CACHE } from '../modules/local/get_db_cache'


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ARETE {

    // Check mandatory parameters
    if (params.input_sample_table) { ch_input = file(params.input_sample_table) } else { exit 1, 'Input samplesheet not specified!' }
    if (params.reference_genome) {
        ch_reference_genome = file(params.reference_genome)
        use_reference_genome = true
        }
    else {
        ch_reference_genome = []
        use_reference_genome = false
    }
    if (!params.skip_poppunk && params.poppunk_model == null) { exit 1, 'A model must be specified with --poppunk_model in order to run PopPunk' }

    ch_bakta_db = params.bakta_db ? file(params.bakta_db) : null
    db_cache = params.db_cache ? params.db_cache : false
    use_full_alignment = params.use_full_alignment
    use_fasttree = params.use_fasttree

    // TODO
    // Outgroup genome isnt currently used for anything. Used to be used for SNIPPY and ended up in the core genome alignment.
    // Look into whether it's possible to include something in the alignment/tree
    //if (params.outgroup_genome ) { ch_outgroup_genome = file(params.outgroup_genome) } else { ch_outgroup_genome = '' }

    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK(ch_input)

    if(db_cache){
        GET_DB_CACHE(db_cache)
        /////////////////// ASSEMBLY ///////////////////////////
        ASSEMBLE_SHORTREADS(
            INPUT_CHECK.out.reads,
            ch_reference_genome,
            use_reference_genome,
            GET_DB_CACHE.out.minikraken
            )
    }
    else{
        ASSEMBLE_SHORTREADS(
            INPUT_CHECK.out.reads,
            ch_reference_genome,
            use_reference_genome,
            []
            )


    }

    ASSEMBLE_SHORTREADS.out.scaffolds.set { assemblies }

    CHECK_ASSEMBLIES (
        assemblies,
        ch_reference_genome,
        use_reference_genome
    )
    ch_software_versions = ch_software_versions.mix(CHECK_ASSEMBLIES.out.assemblyqc_software)

    if (params.apply_filtering) {
        CHECK_ASSEMBLIES.out.assemblies
            .set { assemblies }
    }

    if (db_cache) {
        /////////////////// ANNOTATION ///////////////////////////
        ANNOTATE_ASSEMBLIES(
            assemblies,
            ch_bakta_db,
            GET_DB_CACHE.out.vfdb,
            GET_DB_CACHE.out.cazydb,
            GET_DB_CACHE.out.bacmet,
            GET_DB_CACHE.out.iceberg,
            GET_DB_CACHE.out.card_json,
            GET_DB_CACHE.out.card_version
        )
    } else {
        /////////////////// ANNOTATION ///////////////////////////
        ANNOTATE_ASSEMBLIES(
            assemblies,
            ch_bakta_db,
            [],
            [],
            [],
            [],
            [],
            []
            )
    }

    ANNOTATE_ASSEMBLIES.out.gff.set { gffs }

    if (!params.skip_poppunk) {
        RUN_POPPUNK(assemblies)
        ch_software_versions = ch_software_versions.mix(RUN_POPPUNK.out.poppunk_version)

        if (params.enable_subsetting) {

            SUBSET_GENOMES(
                assemblies,
                RUN_POPPUNK.out.poppunk_db
            )

            SUBSET_GENOMES.out.genomes_to_remove.set { to_remove }

            gffs
                .combine (to_remove.collect().map { [it] })
                .filter { meta, path, to_remove -> !(meta.id in to_remove) }
                .map { it[0, 1] }
                .set { gffs }
        }

        if (params.run_recombination) {
            RECOMBINATION (
                assemblies,
                RUN_POPPUNK.out.clusters,
                CHECK_ASSEMBLIES.out.quast_report
            )
        }
    }

    ch_software_versions = ch_software_versions.mix(ANNOTATE_ASSEMBLIES.out.annotation_software)
    ch_software_versions = ch_software_versions.mix(ASSEMBLE_SHORTREADS.out.assembly_software)

    ////////////////////////// PANGENOME /////////////////////////////////////
    if (!params.skip_phylo) {
        PHYLOGENOMICS(gffs, use_full_alignment, use_fasttree)
        ch_software_versions = ch_software_versions.mix(PHYLOGENOMICS.out.phylo_software)

        if (params.run_evolccm) {
            EVOLCCM (
                PHYLOGENOMICS.out.core_tree,
                ANNOTATE_ASSEMBLIES.out.feature_profile
            )
        }

        if (params.run_rspr) {
            PHYLOGENOMICS.out.gene_trees
                .flatten()
                .map{it -> it.toString() }
                .collectFile(newLine: true) { item ->
                    ["${item}.txt",
                    "sample,path\n" + item + ',' + item ]
                }
                .set { individual_sheets }

            individual_sheets
                .collectFile(name: 'gene_tree_paths.txt', skip:1 , keepHeader: true)
                .set{ gene_tree_sheet }

            RSPR (
                PHYLOGENOMICS.out.core_tree,
                gene_tree_sheet,
                ANNOTATE_ASSEMBLIES.out.annotation
            )
        }
    }

    ////////////////////////// GENE ORDER /////////////////////////////////////
    if (params.run_gene_order) {
        GENE_ORDER (
            ANNOTATE_ASSEMBLIES.out.faa,
            ANNOTATE_ASSEMBLIES.out.gbk,
            ANNOTATE_ASSEMBLIES.out.rgi,
        )
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_software_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    /*
     * MODULE: MultiQC
     */
    workflow_summary    = WorkflowArete.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    //Mix QUAST results into one report file

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(ASSEMBLE_SHORTREADS.out.multiqc)
    ch_multiqc_files = ch_multiqc_files.mix(CHECK_ASSEMBLIES.out.multiqc)
    ch_multiqc_files = ch_multiqc_files.mix(ANNOTATE_ASSEMBLIES.out.multiqc)

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.versions.ifEmpty(null))

}

workflow ASSEMBLY {

    // Check mandatory parameters
    if (params.input_sample_table) { ch_input = file(params.input_sample_table) } else { exit 1, 'Input samplesheet not specified!' }
    if (params.reference_genome) {
        ch_reference_genome = file(params.reference_genome)
        use_reference_genome = true
        }
    else {
        ch_reference_genome = []
        use_reference_genome = false
    }
    //if (params.outgroup_genome ) { ch_outgroup_genome = file(params.outgroup_genome) } else { ch_outgroup_genome = '' }

    // Setup dbcache
    db_cache = params.db_cache ? params.db_cache : false

    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK(ch_input)
    if(db_cache){
        GET_DB_CACHE(db_cache)
        ASSEMBLE_SHORTREADS(INPUT_CHECK.out.reads, ch_reference_genome, use_reference_genome, GET_DB_CACHE.out.minikraken)
    } else {
        ASSEMBLE_SHORTREADS(INPUT_CHECK.out.reads, ch_reference_genome, use_reference_genome, [])
    }

    ch_software_versions = ch_software_versions.mix(ASSEMBLE_SHORTREADS.out.assembly_software)

    CHECK_ASSEMBLIES(
        ASSEMBLE_SHORTREADS.out.scaffolds,
        ch_reference_genome,
        use_reference_genome
    )
    ch_software_versions = ch_software_versions.mix(CHECK_ASSEMBLIES.out.assemblyqc_software)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_software_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //multiqc
    workflow_summary    = WorkflowArete.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(CHECK_ASSEMBLIES.out.multiqc)
    ch_multiqc_files = ch_multiqc_files.mix(ASSEMBLE_SHORTREADS.out.multiqc)

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
}

// annotate existing assemblies
workflow ANNOTATION {
    // Check mandatory parameters
    if (params.input_sample_table) { ch_input = file(params.input_sample_table) } else { exit 1, 'Input samplesheet not specified!' }
    if (params.reference_genome) {
        ch_reference_genome = file(params.reference_genome)
        use_reference_genome = true
        }
    else {
        ch_reference_genome = []
        use_reference_genome = false
    }
    if (!params.skip_poppunk && params.poppunk_model == null) { exit 1, 'A model must be specified with --poppunk_model in order to run PopPunk' }

    ch_bakta_db = params.bakta_db ? file(params.bakta_db) : null

    db_cache = params.db_cache ? params.db_cache : false
    use_full_alignment = params.use_full_alignment ? true : false
    use_fasttree = params.use_fasttree ? true: false

    // TODO
    // Outgroup genome isnt currently used for anything. Used to be used for SNIPPY and ended up in the core genome alignment.
    // Look into whether it's possible to include something in the alignment/tree
    //if (params.outgroup_genome ) { ch_outgroup_genome = file(params.outgroup_genome) } else { ch_outgroup_genome = '' }

    ch_software_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    ANNOTATION_INPUT_CHECK(ch_input)

    ANNOTATION_INPUT_CHECK.out.genomes.set { assemblies }

    if (params.run_recombination || params.apply_filtering) {
        CHECK_ASSEMBLIES(
            assemblies,
            ch_reference_genome,
            use_reference_genome
        )
        ch_software_versions = ch_software_versions.mix(CHECK_ASSEMBLIES.out.assemblyqc_software)
        ch_multiqc_files = ch_multiqc_files.mix(CHECK_ASSEMBLIES.out.multiqc)

        if (params.apply_filtering) {
            CHECK_ASSEMBLIES.out.assemblies
                .set { assemblies }
        }
    }


    if(db_cache){
        GET_DB_CACHE(db_cache)

        /////////////////// ANNOTATION ///////////////////////////
        ANNOTATE_ASSEMBLIES(
            assemblies,
            ch_bakta_db,
            GET_DB_CACHE.out.vfdb,
            GET_DB_CACHE.out.cazydb,
            GET_DB_CACHE.out.bacmet,
            GET_DB_CACHE.out.iceberg,
            GET_DB_CACHE.out.card_json,
            GET_DB_CACHE.out.card_version
            )
    }
    else{

        /////////////////// ANNOTATION ///////////////////////////
        ANNOTATE_ASSEMBLIES(
            assemblies,
            ch_bakta_db,
            [],
            [],
            [],
            [],
            [],
            []
            )
    }
    ch_software_versions = ch_software_versions.mix(ANNOTATE_ASSEMBLIES.out.annotation_software)

    ANNOTATE_ASSEMBLIES.out.gff.set { gffs }

    if (!params.skip_poppunk) {
        RUN_POPPUNK(assemblies)
        ch_software_versions = ch_software_versions.mix(RUN_POPPUNK.out.poppunk_version)

        if (params.enable_subsetting) {

            SUBSET_GENOMES(
                assemblies,
                RUN_POPPUNK.out.poppunk_db
            )

            SUBSET_GENOMES.out.genomes_to_remove.set { to_remove }

            // Filter GFFs not in the to_remove ID list
            gffs
                .combine (to_remove.collect().map { [it] })
                .filter { meta, path, to_remove -> !(meta.id in to_remove) }
                .map { it[0, 1] }
                .set { gffs }
        }

        if (params.run_recombination) {
            RECOMBINATION (
                assemblies,
                RUN_POPPUNK.out.clusters,
                CHECK_ASSEMBLIES.out.quast_report
            )
        }
    }

    ////////////////////////// PANGENOME /////////////////////////////////////
    if (!params.skip_phylo) {
        PHYLOGENOMICS(gffs, use_full_alignment, use_fasttree)
        ch_software_versions = ch_software_versions.mix(PHYLOGENOMICS.out.phylo_software)

        if (params.run_evolccm) {
            EVOLCCM (
                PHYLOGENOMICS.out.core_tree,
                ANNOTATE_ASSEMBLIES.out.feature_profile
            )
        }

        if (params.run_rspr) {
            PHYLOGENOMICS.out.gene_trees
                .flatten()
                .map{it -> it.toString() }
                .collectFile(newLine: true) { item ->
                    ["${item}.txt",
                    "sample,path\n" + item + ',' + item ]
                }
                .set { individual_sheets }

            individual_sheets
                .collectFile(name: 'gene_tree_paths.txt', skip:1 , keepHeader: true)
                .set{ gene_tree_sheet }

            RSPR (
                PHYLOGENOMICS.out.core_tree,
                gene_tree_sheet,
                ANNOTATE_ASSEMBLIES.out.annotation
            )
        }
    }

    ////////////////////////// GENE ORDER /////////////////////////////////////
    if (params.run_gene_order) {
        GENE_ORDER (
            ANNOTATE_ASSEMBLIES.out.faa,
            ANNOTATE_ASSEMBLIES.out.gbk,
            ANNOTATE_ASSEMBLIES.out.rgi,
        )
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_software_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    /*
     * MODULE: MultiQC
     */
    workflow_summary    = WorkflowArete.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    //Mix QUAST results into one report file

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(ANNOTATE_ASSEMBLIES.out.multiqc)

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.versions.ifEmpty(null))

}

workflow QUALITYCHECK {
    if (params.input_sample_table){ ch_input = file(params.input_sample_table) } else { exit 1, 'Input samplesheet not specified!' }
    if (params.reference_genome) {
        ch_reference_genome = file(params.reference_genome)
        use_reference_genome = true
        }
    else {
        ch_reference_genome = []
        use_reference_genome = false
    }
    ch_software_versions = Channel.empty()
    db_cache = params.db_cache ? params.db_cache : false
    ch_multiqc_files = Channel.empty()

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    ANNOTATION_INPUT_CHECK(ch_input)

    CHECK_ASSEMBLIES(
        ANNOTATION_INPUT_CHECK.out.genomes,
        ch_reference_genome,
        use_reference_genome
    )
    ch_software_versions = ch_software_versions.mix(CHECK_ASSEMBLIES.out.assemblyqc_software)
    ch_multiqc_files = ch_multiqc_files.mix(CHECK_ASSEMBLIES.out.multiqc)

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.versions.ifEmpty(null))

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_software_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //multiqc
    workflow_summary    = WorkflowArete.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)
}

workflow PHYLO {
    // Check mandatory parameters
    if (params.input_sample_table) { ch_input = file(params.input_sample_table) } else { exit 1, 'Input samplesheet not specified!' }

    use_full_alignment = params.use_full_alignment ? true : false
    use_fasttree = params.use_fasttree ? true: false

    ch_software_versions = Channel.empty()

    PHYLO_INPUT_CHECK(ch_input)

    PHYLO_INPUT_CHECK.out.genomes.set { gffs }

    ////////////////////////// PANGENOME /////////////////////////////////////
    PHYLOGENOMICS(gffs, use_full_alignment, use_fasttree)
    ch_software_versions = ch_software_versions.mix(PHYLOGENOMICS.out.phylo_software)


    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_software_versions.unique().collectFile(name: 'collated_versions.yml')
    )

}

workflow POPPUNK {

    if (params.input_sample_table) { ch_input = file(params.input_sample_table) } else { exit 1, 'Input samplesheet not specified!' }
    if (!params.skip_poppunk && params.poppunk_model == null) { exit 1, 'A model must be specified with --poppunk_model in order to run PopPunk' }

    /*
    * SUBWORKFLOW: Read in samplesheet, validate and stage input files
    */
    ANNOTATION_INPUT_CHECK(ch_input)

    ANNOTATION_INPUT_CHECK.out.genomes.set { assemblies }

    RUN_POPPUNK(assemblies)
    ch_software_versions = RUN_POPPUNK.out.poppunk_version

    if (params.enable_subsetting) {

        SUBSET_GENOMES(
            assemblies,
            RUN_POPPUNK.out.poppunk_db
        )

    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_software_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}


workflow RUN_RSPR {
    if (params.input_sample_table) { ch_input = Channel.of(file(params.input_sample_table)) } else { exit 1, 'Input samplesheet not specified!' }
    if (params.core_gene_tree) { ch_core = file(params.core_gene_tree) } else { exit 1, 'Core tree not specified!' }
    ch_annotation_data = params.concatenated_annotation ? file(params.concatenated_annotation) : []

    RSPR_INPUT_CHECK (
        ch_input
    )

    RSPR (
        ch_core,
        RSPR_INPUT_CHECK.out.trees,
        ch_annotation_data
    )
}

workflow RUN_EVOLCCM {
    if (params.core_gene_tree) { ch_core = file(params.core_gene_tree) } else { exit 1, 'Core tree not specified!' }
    if (params.feature_profile) { ch_input = file(params.feature_profile) } else { exit 1, 'Input feature profile not specified!' }

    EVOLCCM (
        ch_core,
        ch_input
    )
}

workflow RUN_RECOMBINATION {
    if (params.input_sample_table){ ch_input = file(params.input_sample_table) } else { exit 1, 'Input samplesheet not specified!' }
    if (params.reference_genome) {
        ch_reference_genome = file(params.reference_genome)
        use_reference_genome = true
        }
    else {
        ch_reference_genome = []
        use_reference_genome = false
    }
    if (params.poppunk_model == null) { exit 1, 'A model must be specified with --poppunk_model in order to run PopPunk' }
    ch_software_versions = Channel.empty()

    ANNOTATION_INPUT_CHECK(ch_input)
    ANNOTATION_INPUT_CHECK.out.genomes.set { assemblies }

    CHECK_ASSEMBLIES(
        assemblies,
        ch_reference_genome,
        use_reference_genome
    )
    ch_software_versions = ch_software_versions.mix(CHECK_ASSEMBLIES.out.assemblyqc_software)

    if (params.apply_filtering) {
        CHECK_ASSEMBLIES.out.assemblies
            .set { assemblies }
    }

    RUN_POPPUNK(assemblies)
    ch_software_versions = ch_software_versions.mix(RUN_POPPUNK.out.poppunk_version)


    RECOMBINATION (
        assemblies,
        RUN_POPPUNK.out.clusters,
        CHECK_ASSEMBLIES.out.quast_report
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_software_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    /*
     * MODULE: MultiQC
     */
    workflow_summary    = WorkflowArete.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    //Mix QUAST results into one report file
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(CHECK_ASSEMBLIES.out.multiqc)

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.versions.ifEmpty(null))

}

workflow RUN_GENE_ORDER {
    if (params.input_sample_table){ ch_input = Channel.of(file(params.input_sample_table)) } else { exit 1, 'Input samplesheet not specified!' }

    GENEORDER_INPUT_CHECK (
        ch_input
    )

    GENEORDER_INPUT_CHECK.out.geneorder_input
        .set { all_inputs }

    all_inputs
        .map { it -> [it[0], it[1]] }
        .set { assemblies }

    all_inputs
        .map { it -> [it[0], it[2]] }
        .set { gbks }

    all_inputs
        .map { it -> [it[0], it[3]] }
        .set { rgis }

    GENE_ORDER (
        assemblies,
        gbks,
        rgis
    )
}
/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
