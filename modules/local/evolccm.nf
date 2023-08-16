process EVOLCCM {
    tag "$core_tree"
    label 'process_high'

    conda "bioconda::bioconductor-tximeta=1.8.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/jvfe/evol_ccm:v1.2':
        'docker.io/jvfe/evol_ccm:v1.2' }"

    input:
    path core_tree
    path feature_table

    output:
    path "EvolCCM_*tsv.gz" , emit: profile
    path "EvolCCM_*pvals"  , emit: pvalues
    path "EvolCCM_*X2"     , emit: x2
    path "EvolCCM_*tre"    , emit: tree

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ParallelEvolCCM.R \\
        --intree \\
        $core_tree \\
        --intable \\
        $feature_table \\
        --cores \\
        $task.cpus
    """

    stub:
    """
    touch EvolCCM_test.tsv
    touch EvolCCM_test.tsv.pvals
    touch EvolCCM_test.tsv.X2
    touch EvolCCM_test.tre
    """
}
