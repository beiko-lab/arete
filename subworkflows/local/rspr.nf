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
            params.max_support_threshold,
            params.min_heatmap_approx_rspr,
            params.max_heatmap_approx_rspr,
        )

        RSPR_EXACT (
            RSPR_APPROX.out.csvs.flatten(),
            RSPR_APPROX.out.rooted_reference_tree.first(),
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

        RSPR_EXACT.out.txt
            .collectFile(
                name: 'cluster_file.txt',
                storeDir: "${params.outdir}/dynamics/rSPR/exact"
            )
            .set { cluster_file }

        RSPR_HEATMAP (
            exact_output,
            cluster_file,
            RSPR_APPROX.out.rooted_reference_tree.first(),
            params.min_heatmap_exact_rspr,
            params.max_heatmap_exact_rspr,
        )

    emit:
        exact_output
}
