include { RSPR_APPROX } from '../../modules/local/rspr/approx.nf'
include { RSPR_EXACT } from '../../modules/local/rspr/exact.nf'

workflow RSPR {

    take:
        core_tree
        gene_trees

    main:

        RSPR_APPROX (
            core_tree,
            gene_trees.collect(),
            params.min_rspr_distance
        )

        RSPR_EXACT (
            RSPR_APPROX.out.csvs.flatten(),
            RSPR_APPROX.out.rooted_gene_trees.first()
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
