def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

def modules = params.modules.clone()

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//

include { FASTQC  } from '../../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include { MULTIQC } from '../../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )
include { FASTQC as TRIM_FASTQC } from '../../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include { FASTP                 } from '../../modules/nf-core/modules/fastp/main'  addParams( options: [:] )
include { UNICYCLER             } from '../../modules/nf-core/modules/unicycler/main'  addParams( options: [:] )
include { KRAKEN2_KRAKEN2 as KRAKEN2_RUN } from '../../modules/nf-core/modules/kraken2/kraken2/main' addParams( options: [:] )

include { QUAST                 } from '../../modules/nf-core/modules/quast/main'  addParams( options: [:] )
include { CHECKM_LINEAGEWF } from '../../modules/nf-core/modules/checkm/lineagewf/main' addParams( options: [:] )

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include { KRAKEN2_DB } from '../../modules/local/get_minikraken'  addParams( options: [:] )

// Usage pattern from nf-core/rnaseq: Empty dummy file for optional inputs
ch_dummy_input = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ASSEMBLE_SHORTREADS{
    take:
        reads
        ch_reference_genome
        use_reference_genome
        krakendb_cache

    main:
    /////////////////// Read Processing /////////////////////////////
        /*
        * MODULE: Run FastQC
        */
        ch_software_versions = Channel.empty()
        FASTQC(reads, "raw_fastqc")
        ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))

        /*
        * MODULE: Trim Reads
        */
        FASTP(reads, false, false)
        ch_software_versions = ch_software_versions.mix(FASTP.out.versions.first().ifEmpty(null))

        /*
        * MODULE: Run FastQC on trimmed reads
        */
        TRIM_FASTQC(FASTP.out.reads, "trim_fastqc")
        ch_software_versions = ch_software_versions.mix(TRIM_FASTQC.out.version.first().ifEmpty(null))

        ///*
        // * MODULE: Run Kraken2
        // */
        if (krakendb_cache){
            krakendb_cache.set { ch_kraken_db }
        }
        else {
            KRAKEN2_DB()
            KRAKEN2_DB.out.minikraken.set { ch_kraken_db }
        }


        KRAKEN2_RUN(FASTP.out.reads, ch_kraken_db)
        ch_software_versions = ch_software_versions.mix(KRAKEN2_RUN.out.versions.first().ifEmpty(null))

        /////////////////// ASSEMBLE /////////////////////////////
        /*
        * MODULE: Assembly
        */

        // unicycler can accept short reads and long reads. For now, shortread only: Pass empty list for optional file args
        ch_unicycler_input = FASTP.out.reads.map { it -> it + [[]]}
        UNICYCLER(ch_unicycler_input)
        ch_software_versions = ch_software_versions.mix(UNICYCLER.out.versions.first().ifEmpty(null))

        // Unicycler outputs not quite right for QUAST. Need to re-arrange
        // pattern adapted from nf-core/bacass
        ch_assembly = Channel.empty()
        ch_assembly = ch_assembly.mix(UNICYCLER.out.scaffolds.dump(tag: 'unicycler'))
        ch_assembly
            .map { meta, fasta -> fasta } //QUAST doesn't take the meta tag
            .collect()
            .set { ch_to_quast }
        /*
        * Module: Evaluate Assembly
        */
        QUAST(ch_to_quast, ch_reference_genome, [], use_reference_genome, false)
        ch_software_versions = ch_software_versions.mix(QUAST.out.versions.first().ifEmpty(null))

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(TRIM_FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_RUN.out.txt.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.tsv.collect())

    emit:
        assemblies = ch_assembly
        scaffolds = UNICYCLER.out.scaffolds
        assembly_software = ch_software_versions
        multiqc = ch_multiqc_files

}
