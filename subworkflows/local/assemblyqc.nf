include { GET_DB_CACHE } from '../../modules/local/get_db_cache'
include { QUAST } from '../../modules/nf-core/quast/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_RUN } from '../../modules/nf-core/kraken2/kraken2/main'
include { CHECKM_LINEAGEWF } from '../../modules/nf-core/checkm/lineagewf/main'

workflow CHECK_ASSEMBLIES {
    take:
        assemblies
        reference_genome
        use_reference_genome

    main:

        ch_multiqc_files = Channel.empty()
        ch_software_versions = Channel.empty()

        fasta_extension = assemblies.map{ id, path -> path.getExtension() }.first()

        /*
        * Module: CheckM Quality Check
        */
        CHECKM_LINEAGEWF(assemblies, fasta_extension, [])
        ch_software_versions = ch_software_versions.mix(CHECKM_LINEAGEWF.out.versions.first().ifEmpty(null))
        /*
        * Module: QUAST quality check
        */
        assemblies
            .map { meta, fasta -> fasta.toString() }
            .collectFile(name:'assemblies.txt', newLine: true)
            .set { quast_input }

        QUAST (
            quast_input,
            reference_genome,
            [],
            use_reference_genome,
            false
        )
        ch_software_versions = ch_software_versions.mix(QUAST.out.versions.ifEmpty(null))
        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.tsv.collect())

    emit:
        quast_report = QUAST.out.transposed_report
        assemblyqc_software = ch_software_versions
        multiqc = ch_multiqc_files
}
