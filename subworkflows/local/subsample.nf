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

// workflow {
//     Channel.of(
//                     [[id:'SRR14022735'], "$baseDir/test/SRR14022735_T1.scaffolds.fa"],
//                     [[id:'SRR14022737'], "$baseDir/test/SRR14022737_T1.scaffolds.fa"],
//                     [[id:'SRR14022754'], "$baseDir/test/SRR14022754_T1.scaffolds.fa"],
//                     [[id:'SRR14022764'], "$baseDir/test/SRR14022764_T1.scaffolds.fa"],
//                 )
//                 .set { sample_channel }

//     SUBSET_GENOMES(sample_channel, "../poppunk_bgmm.dists.out")

//     SUBSET_GENOMES.out.filtered_genomes.view()
// }
