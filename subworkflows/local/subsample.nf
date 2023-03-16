workflow SUBSET_GENOMES {

    take:
        genome_assemblies
        poppunk_distances

    main:

        def core_threshold = 100.0 - params.core_similarity
        def accessory_threshold = 100.0 - params.accessory_similarity

        Channel
            .fromPath(poppunk_distances)
            .splitCsv(header: true, sep: '\t')
            .filter { row -> row.Core.toFloat() < core_threshold  && row.Accessory.toFloat() < accessory_threshold}
            .map { row -> row.Query }
            .set { genomes_to_remove }

        genome_assemblies
            .combine (genomes_to_remove.collect().map { [it] })
            .filter { meta, path, to_remove -> !(meta.id in to_remove) }
            .map { it[0, 1] }
            .set { filtered_genomes }

    emit:
        filtered_genomes = filtered_genomes
}
