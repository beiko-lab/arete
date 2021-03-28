// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process QUAST {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id + "/assembly/" + getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? 'bioconda::quast=5.0.2' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl526hb5aa323_2'
    } else {
        container 'quay.io/biocontainers/quast:5.0.2--py37pl526hb5aa323_2'
    }

    input:
    tuple val(meta), path(assembly)
    path(reference_genome)

    output:
    path "${prefix}", emit: results
    path '*.tsv', emit: tsv
    path '*.version.txt', emit: version

    script:
    def software  = getSoftwareName(task.process)
    prefix        = options.suffix ?: software
    """
    quast.py \\
        --output-dir ${prefix} \\
        --threads $task.cpus \\
        -r $reference_genome \\
        $options.args \\
        $assembly
    cp ${prefix}/report.tsv ${meta.id}_report.tsv
    echo \$(quast.py --version 2>&1) | sed 's/^.*QUAST v//; s/ .*\$//' > ${software}.version.txt
    """
}
