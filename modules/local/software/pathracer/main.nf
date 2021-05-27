// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PATHRACER {
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::pathracer=3.15.0.dev" : null)

    input:
    tuple val(meta), path(graph)
    path  hmm

    output:
    tuple val(meta), path('*.all.edges.fa'), emit: all_edges
    tuple val(meta), path('*.log')         , emit: log

    script:
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    pathracer $hmm $graph --output ./ --rescore -t $task.cpus

    mv pathracer.log ${prefix}.pathracer.log
    mv all.edges.fa ${prefix}.all.edges.fa
    """
}
