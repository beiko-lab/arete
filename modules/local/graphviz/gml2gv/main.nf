process GML2GV {
    label 'process_low'

    conda "bioconda::perl-graphviz=2.24"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/perl-graphviz:2.24--pl5321h4b32bfc_1':
        'quay.io/biocontainers/perl-graphviz:2.24--pl5321h4b32bfc_1' }"

    input:
    path(panaroo_results)

    output:
    path("final_graph.dot"), emit: dot
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    gml2gv \\
        $panaroo_results/final_graph.gml \\
        $args \\
        > final_graph.dot

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphviz: \$(echo \$(dot -V 2>&1) | sed 's/dot - graphviz version //')
    END_VERSIONS
    """
    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    touch final_graph.dot

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        graphviz: \$(echo \$(dot -V 2>&1) | sed 's/dot - graphviz version //')
    END_VERSIONS
    """
}
