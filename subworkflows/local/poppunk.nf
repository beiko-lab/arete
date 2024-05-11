include { POPPUNK_MAKE_SAMPLESHEET } from '../../modules/local/poppunk_samplesheet'
include { POPPUNK_CREATEDB } from '../../modules/local/poppunk/createdb/main'
include { POPPUNK_QCDB } from '../../modules/local/poppunk/qcdb/main'
include { POPPUNK_FITMODEL } from '../../modules/local/poppunk/fitmodel/main'
include { POPPUNK_VISUALISE } from '../../modules/local/poppunk/visualise/main'

workflow RUN_POPPUNK {

    take:
        genome_assemblies

    main:

        if (params.reference_genome && params.run_poppunk_qc) {
            Channel.fromList([
                [[id:'reference_genome'], file(params.reference_genome)]
            ])
                .set { ref_genome }

            genome_assemblies
                .mix(ref_genome)
                .map { it -> [ it[0], it[1] ] }
                .set { genome_assemblies }
        }

        genome_assemblies
            .map { meta, path -> [meta.id, path.getName()] }
            .collectFile(newLine: true) { item ->
                [ "${item[0]}.txt", item[0] + '\t' + item[1] ]
            }
            .collectFile(name: 'poppunk_samplesheet.tsv')
            .set { poppunk_samplesheet }

        POPPUNK_CREATEDB(poppunk_samplesheet, genome_assemblies.collect { it[1] } )

        POPPUNK_CREATEDB.out.poppunk_db.set { poppunk_db }

        if (params.run_poppunk_qc) {
            type_isolate = params.reference_genome ? "reference_genome" : []
            POPPUNK_QCDB(poppunk_db, type_isolate)
            POPPUNK_QCDB.out.poppunk_db.set{ poppunk_db }
        }

        POPPUNK_FITMODEL(poppunk_db, params.poppunk_model)

        POPPUNK_FITMODEL.out.poppunk_results.set { poppunk_results }

        POPPUNK_VISUALISE(poppunk_results, poppunk_db)

    emit:
        poppunk_version = POPPUNK_FITMODEL.out.versions.ifEmpty(null)
        poppunk_results = poppunk_results
        poppunk_visualisations = POPPUNK_VISUALISE.out.poppunk_visualizations
        clusters = POPPUNK_FITMODEL.out.clusters
        poppunk_db = poppunk_db
}
