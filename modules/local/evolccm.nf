process EVOLCCM {
    tag "$core_tree"
    label 'process_high'
    label 'process_long'

    conda "bioconda::bioconductor-tximeta=1.8.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/jvfe/evol_ccm:v1.3.1':
        'docker.io/jvfe/evol_ccm:v1.3.1' }"

    input:
    path core_tree
    path feature_table

    output:
    path "EvolCCM_*tsv.gz"     , emit: profile
    path "EvolCCM_*_pvals.tsv" , emit: pvalues
    path "EvolCCM_*_X2.tsv"    , emit: x2
    path "EvolCCM_*tre"        , emit: tree

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ParallelEvolCCM.R \\
        --intree \\
        $core_tree \\
        --intable \\
        $feature_table \\
        --cores \\
        $task.cpus \\
        $args
    """

    stub:
    """
    touch EvolCCM_test.tsv.gz
    touch EvolCCM_test_pvals.tsv
    touch EvolCCM_test_X2.tsv
    touch EvolCCM_test.tre
    """
}
