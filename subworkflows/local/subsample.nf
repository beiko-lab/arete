include { POPPUNK_EXTRACT_DISTANCES } from '../../modules/local/poppunk/extractdistances/main'
include { MAKE_HEATMAP_AND_FILTER } from "../../modules/local/make_heatmap_and_filter.nf"

workflow SUBSET_GENOMES {

    take:
        genome_assemblies
        poppunk_db

    main:
        def accessory_threshold = (100.0 - params.accessory_similarity)/100.0
        def core_threshold = (100.0 - params.core_similarity)/100.0

        println '\033[0;33mWARN: Your assemblies will be subsampled (--enable_subsetting is true)\033[0m'

        POPPUNK_EXTRACT_DISTANCES(poppunk_db)
        POPPUNK_EXTRACT_DISTANCES.out.poppunk_distances.set{ poppunk_distances }

        MAKE_HEATMAP_AND_FILTER (
            poppunk_distances,
            core_threshold,
            accessory_threshold
        )

        MAKE_HEATMAP_AND_FILTER.out.removed_genomes
            .splitCsv(sep: '\t')
            .map { it }
            .set { genomes_to_remove }

    emit:
        genomes_to_remove = genomes_to_remove
}
