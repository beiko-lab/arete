include { FILTER_AND_CONCAT_MATCHES } from '../../modules/local/filter_matches'


workflow FILTER_ALIGNMENT {

    take:
        diamond_results
        db_name
        blast_columns
        min_pident
        min_qcover

    main:
        diamond_results
            .filter{ id, path -> path.size() > 0 }
            .collect { id, path -> path }
            .set { results }

        FILTER_AND_CONCAT_MATCHES(
            results,
            db_name,
            blast_columns,
            min_pident,
            min_qcover,
            1,
            false
        )
        FILTER_AND_CONCAT_MATCHES.out.txt.set { diamond_filtered }

    emit:
        concatenated = diamond_filtered
}
