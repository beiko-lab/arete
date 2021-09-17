// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ROARY {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::roary=3.13.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/roary:3.13.0--pl526h516909a_0"
    } else {
        container "quay.io/biocontainers/roary:3.13.0--pl526h516909a_0"
    }

    input:
    // TODO check input
    tuple val(meta), path(gff)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("${prefix}/*.bam"), emit: 
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    
    //TODO manage output?
    """
    roary \\
        $options.args \\
        -p $task.cpus \\ //n threads
        -f ${prefix} \\ // -f specifies output directory
        $gff

    echo \$(roary -w 2>&1) | sed 's/^.*roary //; s/Using.*\$//' > ${software}.version.txt
    """
}
