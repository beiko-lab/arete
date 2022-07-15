// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process POPPUNK_MAKESHEET{

}

process POPPUNK_MAKEDB {
    tag "$meta.id"
    label 'process_medium'
    
    conda (params.enable_conda ? "bioconda::poppunk=2.4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.4.0--py39h8884e85_2':
        'quay.io/biocontainers/poppunk:2.4.0--py38hf6d6cf9_1' }"

    input:

    tuple val(meta), path(filesheet)
    val(outdir)

    output:
    
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    samtools \\
        sort \\
        $args \\
        -@ $task.cpus \\
        -o ${prefix}.bam \\
        -T $prefix \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        poppunk: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}

process POPPUNK_RUN{

}