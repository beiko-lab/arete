include { POPPUNK_EXTRACT_DISTANCES } from '../../modules/local/poppunk/extractdistances/main'

workflow SUBSET_GENOMES {

    take:
        genome_assemblies
        poppunk_db

    main:

        POPPUNK_EXTRACT_DISTANCES(poppunk_db)
        POPPUNK_EXTRACT_DISTANCES.out.poppunk_distances.set{ poppunk_distances }

        def core_threshold = 100.0 - params.core_similarity
        def accessory_threshold = 100.0 - params.accessory_similarity

        poppunk_distances
            .splitCsv(header: true, sep: '\t')
            .filter { row -> (row.Core.toFloat() * 100) < core_threshold  && (row.Accessory.toFloat() * 100) < accessory_threshold}
            .map { row -> row.Query }
            .set { genomes_to_remove }

        genomes_to_remove
            .unique()
            .collectFile(newLine: true)
            .collectFile(name: 'removed_genomes.txt', storeDir: "${params.outdir}/poppunk_results")

        genome_assemblies
            .combine (genomes_to_remove.collect().map { [it] })
            .filter { meta, path, to_remove -> !(meta.id in to_remove) }
            .map { it[0, 1] }
            .set { filtered_genomes }

    emit:
        filtered_genomes = filtered_genomes
}
