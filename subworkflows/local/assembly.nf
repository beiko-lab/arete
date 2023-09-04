//
// MODULE: Installed directly from nf-core/modules
//

include { FASTQC as RAW_FASTQC  } from '../../modules/nf-core/fastqc/main'
include { FASTQC as TRIM_FASTQC } from '../../modules/nf-core/fastqc/main'
include { FASTP                 } from '../../modules/nf-core/fastp/main'
include { UNICYCLER             } from '../../modules/nf-core/unicycler/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_RUN } from '../../modules/nf-core/kraken2/kraken2/main'

include { QUAST                 } from '../../modules/nf-core/quast/main'
include { CHECKM_LINEAGEWF } from '../../modules/nf-core/checkm/lineagewf/main'

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include { KRAKEN2_DB } from '../../modules/local/get_minikraken'  addParams( options: [:] )

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
        ch_multiqc_files = Channel.empty()
    /////////////////// Read Processing /////////////////////////////
        /*
        * MODULE: Run FastQC
        */
        ch_software_versions = Channel.empty()
        RAW_FASTQC(reads)
        ch_software_versions = ch_software_versions.mix(RAW_FASTQC.out.versions.first().ifEmpty(null))

        /*
        * MODULE: Trim Reads
        */
        FASTP(reads, [], false, false)
        ch_software_versions = ch_software_versions.mix(FASTP.out.versions.first().ifEmpty(null))

        /*
        * MODULE: Run FastQC on trimmed reads
        */
        TRIM_FASTQC(FASTP.out.reads)
        ch_software_versions = ch_software_versions.mix(TRIM_FASTQC.out.versions.first().ifEmpty(null))

        ///*
        // * MODULE: Run Kraken2
        // */
        if (!params.skip_kraken) {
            if (krakendb_cache){
                krakendb_cache.set { ch_kraken_db }
            }
            else {
                KRAKEN2_DB()
                KRAKEN2_DB.out.minikraken.set { ch_kraken_db }
            }


            KRAKEN2_RUN(FASTP.out.reads, ch_kraken_db, true, true)
            ch_software_versions = ch_software_versions.mix(KRAKEN2_RUN.out.versions.first().ifEmpty(null))
            ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_RUN.out.report.collect{it[1]}.ifEmpty([]))
        }
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
        ch_multiqc_files = ch_multiqc_files.mix(RAW_FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(TRIM_FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    emit:
        scaffolds = UNICYCLER.out.scaffolds
        assembly_software = ch_software_versions
        multiqc = ch_multiqc_files

}
