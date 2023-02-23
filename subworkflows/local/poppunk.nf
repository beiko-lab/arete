include { POPPUNK_MAKE_SAMPLESHEET } from '../../modules/local/poppunk'
include { POPPUNK_MAKEDB } from '../../modules/local/poppunk'
include { POPPUNK_FITMODEL } from '../../modules/local/poppunk'

workflow RUN_POPPUNK {
    Channel
        .fromPath("../assemblies/ann_samplesheet.csv")
        .splitCsv(header: true)
        .map { row -> [[id: row.sample], file(row.fna_file_path)]}
        .set { scaffolds }

    scaffolds
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

    POPPUNK_FITMODEL(poppunk_db)
}
