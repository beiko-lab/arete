include { GET_DB_CACHE } from '../../modules/local/get_db_cache'
include { QUAST } from '../../modules/nf-core/quast/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_RUN } from '../../modules/nf-core/kraken2/kraken2/main'
include { CHECKM_LINEAGEWF } from '../../modules/nf-core/checkm/lineagewf/main'

workflow CHECK_ASSEMBLIES {
    take:
        assemblies
        krakendb_cache
        reference_genome
        use_reference_genome

    main:

        ch_multiqc_files = Channel.empty()
        ch_software_versions = Channel.empty()

        ///*
        // * MODULE: Run Kraken2
        // */
        if (!params.skip_kraken) {
            if(krakendb_cache) {
                GET_DB_CACHE(krakendb_cache)
                KRAKEN2_RUN(assemblies, GET_DB_CACHE.out.minikraken, false, true)
            } else {
                KRAKEN2_DB()
                KRAKEN2_RUN(assemblies, KRAKEN2_DB.out.minikraken, false, true)
            }

            ch_software_versions = ch_software_versions.mix(KRAKEN2_RUN.out.versions.first().ifEmpty(null))
            ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_RUN.out.report.collect{it[1]}.ifEmpty([]))
        }

        fasta_extension = assemblies.map{ id, path -> path.getExtension() }.first()

        /*
        * Module: CheckM Quality Check
        */
        CHECKM_LINEAGEWF(assemblies, fasta_extension, [])
        ch_software_versions = ch_software_versions.mix(CHECKM_LINEAGEWF.out.versions.first().ifEmpty(null))
        /*
        * Module: QUAST quality check
        */
        // Need to reformat assembly channel for QUAST
        // pattern adapted from nf-core/bacass
        ch_assembly = Channel.empty()
        ch_assembly = ch_assembly.mix(assemblies.dump(tag: 'assembly'))
        ch_assembly
            .map { meta, fasta -> fasta } //QUAST doesn't take the meta tag
            .collect()
            .set { ch_to_quast }
        QUAST(ch_to_quast, reference_genome, [], use_reference_genome, false)
        ch_software_versions = ch_software_versions.mix(QUAST.out.versions.ifEmpty(null))

        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.tsv.collect())

    emit:
        assemblyqc_software = ch_software_versions
        multiqc = ch_multiqc_files
}
