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

        RSPR_EXACT (
            RSPR_APPROX.out.csvs.flatten(),
            RSPR_APPROX.out.rooted_gene_trees.first(),
            params.min_branch_length,
            params.max_support_threshold,
            params.max_approx_rspr
        )

        RSPR_EXACT.out.tsv
            .collectFile(
                name: 'exact_output.tsv',
                keepHeader: true,
                storeDir: "${params.outdir}/dynamics/rSPR/exact",
                skip: 1
            )
            .set { exact_output }

        RSPR_HEATMAP (
            exact_output
        )

    emit:
        exact_output
}
