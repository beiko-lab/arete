
/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// TODO nf-core: Add all file path parameters for the pipeline to the lis below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input_sample_table, params.multiqc_config, params.reference_genome ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input_sample_table) { ch_input = file(params.input_sample_table) } else { exit 1, 'Input samplesheet not specified!' }
if (!params.skip_poppunk && params.poppunk_model == null) { exit 1, 'A model must be specified with --poppunk_model in order to run PopPunk' }


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
include { INPUT_CHECK;
          ANNOTATION_INPUT_CHECK } from '../subworkflows/local/input_check'

include { ASSEMBLE_SHORTREADS } from '../subworkflows/local/assembly'
include { ANNOTATE_ASSEMBLIES } from '../subworkflows/local/annotation'
include { PHYLOGENOMICS } from '../subworkflows/local/phylo'
include { RUN_POPPUNK } from '../subworkflows/local/poppunk'
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
include { QUAST                 } from '../modules/nf-core/quast/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_RUN } from '../modules/nf-core/kraken2/kraken2/main'
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
include { CHECKM_LINEAGEWF } from '../modules/nf-core/checkm/lineagewf/main'
//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include { RGI;
          UPDATE_RGI_DB } from '../modules/local/rgi'
include { MOB_RECON } from '../modules/local/mobsuite'
include { KRAKEN2_DB } from '../modules/local/get_minikraken'
include { GET_DB_CACHE } from '../modules/local/get_db_cache'


// Usage pattern from nf-core/rnaseq: Empty dummy file for optional inputs
ch_dummy_input = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

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

        /////////////////// ANNOTATION ///////////////////////////
        ANNOTATE_ASSEMBLIES(
            ASSEMBLE_SHORTREADS.out.scaffolds,
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
        ASSEMBLE_SHORTREADS(
            INPUT_CHECK.out.reads,
            ch_reference_genome,
            use_reference_genome,
            []
            )

        /////////////////// ANNOTATION ///////////////////////////
        ANNOTATE_ASSEMBLIES(
            ASSEMBLE_SHORTREADS.out.scaffolds,
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
    ch_software_versions = ch_software_versions.mix(ASSEMBLE_SHORTREADS.out.assembly_software)

    if (!params.skip_poppunk) {
        RUN_POPPUNK(ASSEMBLE_SHORTREADS.out.scaffolds)
        ch_software_versions = ch_software_versions.mix(RUN_POPPUNK.out.poppunk_version)
    }

    ////////////////////////// PANGENOME /////////////////////////////////////
    PHYLOGENOMICS(ANNOTATE_ASSEMBLIES.out.gff, use_full_alignment, use_fasttree)
    ch_software_versions = ch_software_versions.mix(PHYLOGENOMICS.out.phylo_software)

    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }
    GET_SOFTWARE_VERSIONS(ch_software_versions)

    /*
     * MODULE: MultiQC
     */
    workflow_summary    = WorkflowArete.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    //Mix QUAST results into one report file

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(ASSEMBLE_SHORTREADS.out.multiqc)
    //ch_multiqc_files = ch_multiqc_files.mix(PROKKA.out.txt.collect{it[1]}.ifEmpty([]))
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

    // Get unique list of files containing version information
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }
    GET_SOFTWARE_VERSIONS(ch_software_versions)

    //multiqc
    workflow_summary    = WorkflowArete.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ASSEMBLE_SHORTREADS.out.multiqc)
    //ch_multiqc_files = ch_multiqc_files.mix(PROKKA.out.txt.collect{it[1]}.ifEmpty([]))

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.versions.ifEmpty(null))

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

    ch_bakta_db = params.bakta_db ? file(params.bakta_db) : null

    //db_cache = params.db_cache ? params.db_cache: false
    //ch_db_cache = Channel.empty()
    // ch_assembly_db_cache = Channel.empty()
    // ch_annotation_db_cache = Channel.empty()
    // if (params.db_cache){
    //     ch_assembly_db_cache = GET_ASSEMBLY_DB_CACHE(file(params.db_cache))
    //     ch_annotation_db_cache = GET_ANNOTATION_DB_CACHE(file(params.db_cache))
    // }
    // else{
    //     ch_assembly_db_cache = false
    //     ch_annotation_db_cache = false
    // }
    db_cache = params.db_cache ? params.db_cache : false
    use_full_alignment = params.use_full_alignment ? true : false
    use_fasttree = params.use_fasttree ? true: false

    // TODO
    // Outgroup genome isnt currently used for anything. Used to be used for SNIPPY and ended up in the core genome alignment.
    // Look into whether it's possible to include something in the alignment/tree
    //if (params.outgroup_genome ) { ch_outgroup_genome = file(params.outgroup_genome) } else { ch_outgroup_genome = '' }

    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    ANNOTATION_INPUT_CHECK(ch_input)

    if(db_cache){
        GET_DB_CACHE(db_cache)

        /////////////////// ANNOTATION ///////////////////////////
        ANNOTATE_ASSEMBLIES(
            ANNOTATION_INPUT_CHECK.out.genomes,
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
            ANNOTATION_INPUT_CHECK.out.genomes,
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

    if (!params.skip_poppunk) {
        RUN_POPPUNK(ANNOTATION_INPUT_CHECK.out.genomes)
        ch_software_versions = ch_software_versions.mix(RUN_POPPUNK.out.poppunk_version)
    }

    ////////////////////////// PANGENOME /////////////////////////////////////
    PHYLOGENOMICS(ANNOTATE_ASSEMBLIES.out.gff, use_full_alignment, use_fasttree)
    ch_software_versions = ch_software_versions.mix(PHYLOGENOMICS.out.phylo_software)

    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }
    GET_SOFTWARE_VERSIONS(ch_software_versions)

    /*
     * MODULE: MultiQC
     */
    workflow_summary    = WorkflowArete.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    //Mix QUAST results into one report file

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    //ch_multiqc_files = ch_multiqc_files.mix(PROKKA.out.txt.collect{it[1]}.ifEmpty([]))
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

workflow QUALITYCHECK{
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

    ///*
    // * MODULE: Run Kraken2
    // */
    if (!params.skip_kraken) {
        if(db_cache) {
            GET_DB_CACHE(db_cache)
            KRAKEN2_RUN(ANNOTATION_INPUT_CHECK.out.genomes, GET_DB_CACHE.out.minikraken, false, true)
        } else {
            KRAKEN2_DB()
            KRAKEN2_RUN(ANNOTATION_INPUT_CHECK.out.genomes, KRAKEN2_DB.out.minikraken, false, true)
        }

        ch_software_versions = ch_software_versions.mix(KRAKEN2_RUN.out.versions.first().ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_RUN.out.report.collect{it[1]}.ifEmpty([]))
    }
    /*
    * Module: CheckM Quality Check
    */
    CHECKM_LINEAGEWF(ANNOTATION_INPUT_CHECK.out.genomes, "fna", []) //todo figure out a way to infer the file extension during input check
    ch_software_versions = ch_software_versions.mix(CHECKM_LINEAGEWF.out.versions.first().ifEmpty(null))
    /*
     * Module: QUAST quality check
     */
    // Need to reformat assembly channel for QUAST
    // pattern adapted from nf-core/bacass
    ch_assembly = Channel.empty()
    ch_assembly = ch_assembly.mix(ANNOTATION_INPUT_CHECK.out.genomes.dump(tag: 'assembly'))
    ch_assembly
        .map { meta, fasta -> fasta } //QUAST doesn't take the meta tag
        .collect()
        .set { ch_to_quast }
    QUAST(ch_to_quast, ch_reference_genome, [], use_reference_genome, false)
    ch_software_versions = ch_software_versions.mix(QUAST.out.versions.ifEmpty(null))

    ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.tsv.collect())

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.collect().ifEmpty([]),
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_multiqc_logo.collect().ifEmpty([])
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.versions.ifEmpty(null))

    // Get unique list of files containing version information
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }
    GET_SOFTWARE_VERSIONS(ch_software_versions)
    //multiqc
    workflow_summary    = WorkflowArete.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)
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
