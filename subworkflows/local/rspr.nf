include { RSPR_APPROX } from '../../modules/local/rspr/approx.nf'
include { RSPR_EXACT } from '../../modules/local/rspr/exact.nf'
include { RSPR_HEATMAP } from '../../modules/local/rspr/heatmap.nf'

workflow RSPR {

    take:
        core_tree
        gene_tree_sheet
        annotation

    main:

        RSPR_APPROX (
            core_tree,
            gene_tree_sheet,
            annotation,
            params.min_rspr_distance,
            params.min_branch_length,
            params.max_support_threshold
        )
}
