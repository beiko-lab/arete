// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_NCBI_AMR_HMM {
    tag "NCBI AMRFinder HMM download"
    
    output:
    path "AMR.LIB", emit: hmm

    script:
    """
    wget https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinder/data/latest/AMR.LIB
    """
}

process PATHRACER {
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id + "/annotation/" + getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::pathracer=3.15.0.dev" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pathracer:3.15.0.dev--h2d02072_0"
    } else {
        container "quay.io/biocontainers/pathracer:3.15.0.dev--h2d02072_0"
    }
 
    input:
    tuple val(meta), path(graph)
    path hmm

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
