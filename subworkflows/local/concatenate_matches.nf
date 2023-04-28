include { CONCAT_OUTPUT } from '../../modules/local/concat_output'
include { FILTER_MATCHES } from '../../modules/local/filter_matches'


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
            .set { results }

        FILTER_MATCHES(
            results,
            db_name,
            blast_columns,
            min_pident,
            min_qcover
        )
        FILTER_MATCHES.out.txt.set { diamond_filtered }

        diamond_filtered
            .collect{ id, paths -> paths }
            .set { paths }

        CONCAT_OUTPUT(paths, db_name, 1)
        CONCAT_OUTPUT.out.txt.set { concatenated }

    emit:
        concatenated = concatenated
}
