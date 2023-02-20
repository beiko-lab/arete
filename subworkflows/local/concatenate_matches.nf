include { ADD_COLUMN } from '../../modules/local/add_column.nf'
include { CONCAT_ALIGNMENT } from '../../modules/local/concat_alignments.nf'


workflow FILTER_ALIGNMENT {

    take:
        diamond_results
        db_name
        blast_columns

    main:
        diamond_results
            .filter{ id, path -> path.size() > 0 }
            .set { results }

        ADD_COLUMN(results, db_name)
        ADD_COLUMN.out.txt.set { diamond_added_column }

        diamond_added_column
            .collect()
            .set { paths }

        CONCAT_ALIGNMENT(paths, db_name, blast_columns)
        CONCAT_ALIGNMENT.out.txt.set { concatenated }

    emit:
        concatenated = concatenated
}
