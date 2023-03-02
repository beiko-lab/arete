include { POPPUNK_MAKE_SAMPLESHEET } from '../../modules/local/poppunk_samplesheet'
include { POPPUNK_MAKEDB } from '../../modules/local/poppunk/makedb/main'
include { POPPUNK_QCDB } from '../../modules/local/poppunk/qcdb/main'
include { POPPUNK_FITMODEL } from '../../modules/local/poppunk/fitmodel/main'
include { POPPUNK_VISUALISE } from '../../modules/local/poppunk/visualise/main'

workflow RUN_POPPUNK {

    take:
        genome_assemblies

    main:

        genome_assemblies
            .map { meta, path -> [meta.id, path.toString()] }
            .collectFile(newLine: true) { item ->
                [ "${item[0]}.txt", item[0] + '\t' + item[1] ]
            }
            .collect()
            .set { samplesheets }

        POPPUNK_MAKE_SAMPLESHEET(samplesheets)

        POPPUNK_MAKE_SAMPLESHEET.out.full_samplesheet.set { poppunk_samplesheet }

        POPPUNK_MAKEDB(poppunk_samplesheet)

        POPPUNK_MAKEDB.out.poppunk_db.set { poppunk_db }

        if (params.run_poppunk_qc) {
            POPPUNK_QCDB(poppunk_db, [])
            POPPUNK_QCDB.out.poppunk_db.set{ poppunk_db }
        }

        POPPUNK_FITMODEL(poppunk_db, params.poppunk_model)

        POPPUNK_FITMODEL.out.poppunk_results.set { poppunk_results }

        POPPUNK_VISUALISE(poppunk_results, poppunk_db)

    emit:
        poppunk_version = POPPUNK_FITMODEL.out.versions.ifEmpty(null)
        poppunk_results = poppunk_results
        poppunk_visualisations = POPPUNK_VISUALISE.out.poppunk_visualizations
}
