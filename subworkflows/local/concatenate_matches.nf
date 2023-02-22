include { ADD_COLUMN } from '../../modules/local/add_column'
include { CONCAT_ALIGNMENT } from '../../modules/local/concat_alignments'
include { FILTER_MATCHES } from '../../modules/local/filter_matches'


workflow FILTER_ALIGNMENT {

    take:
        diamond_results
        db_name
        blast_columns

    main:
        diamond_results
            .filter{ id, path -> path.size() > 0 }
            .set { results }

        FILTER_MATCHES(results, db_name, blast_columns)
        FILTER_MATCHES.out.txt.set { diamond_filtered }

        diamond_filtered
            .collect{ id, paths -> paths }
            .set { paths }

        CONCAT_ALIGNMENT(paths, db_name)
        CONCAT_ALIGNMENT.out.txt.set { concatenated }

    emit:
        concatenated = concatenated
}
