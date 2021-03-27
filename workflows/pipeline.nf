////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Validate input parameters
Workflow.validateWorkflowParams(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input_sample_table, params.multiqc_config, params.reference_genome ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input_sample_table) { ch_input = file(params.input_sample_table) } else { exit 1, 'Input samplesheet not specified!' }

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --       IMPORT MODULES / SUBWORKFLOWS      -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''

// Modules: local
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions'   addParams( options: [publish_files : ['csv':'']] )

// Modules: nf-core/modules
include { FASTQC                } from '../modules/nf-core/software/fastqc/main'  addParams( options: modules['fastqc']            )
include { FASTQC as TRIM_FASTQC } from '../modules/nf-core/software/fastqc/main'  addParams( options: modules['fastqc']            )
include { FASTP                 } from '../modules/nf-core/software/fastp/main'  addParams( options: [:]                          )
include { UNICYCLER             } from '../modules/nf-core/software/unicycler/main'  addParams( options: [:] )
include { QUAST                 } from '../modules/nf-core/software/quast/main'  addParams( options: [:]                          )
include { MULTIQC               } from '../modules/nf-core/software/multiqc/main' addParams( options: multiqc_options              )

// Subworkflows: local
include { INPUT_CHECK           } from '../subworkflows/local/input_check'        addParams( options: [:]                          )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report = []

workflow ARETE {

    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    INPUT_CHECK(ch_input)

    /*
     * MODULE: Run FastQC
     */
    FASTQC(INPUT_CHECK.out.reads)
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))
    
    /*
     * MODULE: Trim Reads
     */
    FASTP(INPUT_CHECK.out.reads)
    ch_software_versions = ch_software_versions.mix(FASTP.out.version.first().ifEmpty(null))

    /*
     * MODULE: Run FastQC on trimmed reads
     */
    TRIM_FASTQC(FASTP.out.reads)
    ch_software_versions = ch_software_versions.mix(TRIM_FASTQC.out.version.first().ifEmpty(null))

    ///*
    // * MODULE: Run Kraken2
    // */
    //KRAKEN2(FASTP.out.reads)
    //ch_software_versions = ch_software_versions.mix(KRAKEN2.out.version.first().ifEmpty(null))

    /*
     * MODULE: Assembly
     */
    UNICYCLER(FASTP.out.reads)
    ch_software_versions = ch_software_versions.mix(UNICYCLER.out.version.first().ifEmpty(null))

    /*
     * Module: Evaluate Assembly
     */
    QUAST(UNICYCLER.out.scaffolds)
    ch_software_versions = ch_software_versions.mix(QUAST.out.version.first().ifEmpty(null))

    ///*
    // * Module: Annotate AMR
    // */
    //RGI(UNICYCLER.out.scaffolds)
    //ch_software_versions = ch_software_versions.mix(RGI.out.version.first().ifEmpty(null))

    ///*
    // * Module: Annotate VF
    // */
    //ABRICATE(UNICYCLER.out.assemblies, "VFDB")
    //ch_software_versions = ch_software_versions.mix(ABRICATE.out.version.first().ifEmpty(null))

    ///*
    // * Module: Annotate BacMet
    // */
    //ABRICATE(UNICYCLER.out.assemblies, "BacMet2")
    //ch_software_versions = ch_software_versions.mix(ABRICATE.out.version.first().ifEmpty(null))

    ///*
    // * Module: Prokka
    // */
    //PROKKA(UNICYCLER.out.assemblies)
    //ch_software_versions = ch_software_versions.mix(PROKKA.out.version.first().ifEmpty(null))

    ///*
    // * Module: Mob-Suite
    // */
    //MOBSUITE(UNICYCLER.out.assemblies)
    //ch_software_versions = ch_software_versions.mix(MOBSUITE.out.version.first().ifEmpty(null))

    ///*
    // * Module: Snippy
    // */
    //SNIPPY(UNICYCLER.out.assemblies, ch_reference_genome)
    //ch_software_versions = ch_software_versions.mix(QUAST.out.version.first().ifEmpty(null))
    //
    ///*
    // * Module: IQTree
    // */
    //IQTREE(SNIPPY.out.snps)
    //ch_software_versions = ch_software_versions.mix(IQTREE.out.version.first().ifEmpty(null))


    /*
     * MODULE: Pipeline reporting
     */
    // Get unique list of files containing version information
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
    workflow_summary    = Workflow.paramsSummaryMultiqc(workflow, params.summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(TRIM_FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    //ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2.out.zip.collect{it[1]}.ifEmpty([]))
    //ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.zip.collect{it[1]}.ifEmpty([]))
    
    MULTIQC (ch_multiqc_files.collect())
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
    
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
///////////////////////////////////////////////////
