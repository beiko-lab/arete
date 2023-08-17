include { EVOLCCM as RUN_EVOLCCM } from '../../modules/local/evolccm.nf'

workflow EVOLCCM {

    take:
        core_tree
        feature_profile

    main:

        RUN_EVOLCCM (
            core_tree,
            feature_profile
        )

    emit:
        profile = RUN_EVOLCCM.out.profile
        pvalues = RUN_EVOLCCM.out.pvalues
        x2 = RUN_EVOLCCM.out.x2
        tree = RUN_EVOLCCM.out.tree
}
