include { RSPR_APPROX } from '../../modules/local/rspr/approx.nf'
include { RSPR_EXACT } from '../../modules/local/rspr/exact.nf'

workflow RSPR {

    take:
        core_tree
        gene_trees

    main:

        gene_trees
            .flatten()
            .map{it -> it.toString() }
            .collectFile(name: 'gene_tree_paths.txt', newLine: true)
            .set{ gene_tree_sheet }

        RSPR_APPROX (
            core_tree,
            gene_tree_sheet,
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

        RSPR_EXACT.out.csv
            .collectFile(
                name: 'exact_output.csv',
                keepHeader: true,
                storeDir: "${params.outdir}/dynamics/rSPR/exact",
                skip: 1
            )
            .set { exact_output }

    emit:
        exact_output
}
