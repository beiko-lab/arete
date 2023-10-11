include { POPPUNK_EXTRACT_DISTANCES } from '../../modules/local/poppunk/extractdistances/main'
include { MAKE_HEATMAP } from "../../modules/local/make_heatmap.nf"

workflow SUBSET_GENOMES {

    take:
        genome_assemblies
        poppunk_db

    main:

        println '\033[0;33mWARN: Your assemblies will be subsampled (--enable_subsetting is true)\033[0m'

        POPPUNK_EXTRACT_DISTANCES(poppunk_db)
        POPPUNK_EXTRACT_DISTANCES.out.poppunk_distances.set{ poppunk_distances }

        MAKE_HEATMAP(poppunk_distances)

        def core_threshold = 100.0 - params.core_similarity
        def accessory_threshold = 100.0 - params.accessory_similarity

        poppunk_distances
            .splitCsv(header: true, sep: '\t')
            .filter { row -> (row.Core.toFloat() * 100) < core_threshold  && (row.Accessory.toFloat() * 100) < accessory_threshold }
            .map { row -> row.Query }
            .set { genomes_to_remove }

        genomes_to_remove
            .unique()
            .collectFile(newLine: true)
            .collectFile(name: 'removed_genomes.txt', storeDir: "${params.outdir}/poppunk_results")

    emit:
        genomes_to_remove = genomes_to_remove
}
