include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process UPDATE_RGI_DB {
    tag "CARD"
    label 'process_low'
    conda (params.enable_conda ? "bioconda::rgi=5.1.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/5.1.1-0"
    } else {
        container "quay.io/biocontainers/rgi:5.1.1--py_0"
    }

    output:
    path "card.version.txt", emit: version

    script:
    """
    curl https://card.mcmaster.ca/latest/data --output card.tar.bz2
    tar xvf card.tar.bz2
    rgi load -i card.json
    rgi database -v > card.version.txt
    """
}

process RGI {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id + "/annotation/" + getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::rgi=5.1.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/5.1.1-0"
    } else {
        container "quay.io/biocontainers/rgi:5.1.1--py_0"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_rgi.txt"), emit: tsv
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    prefix = "${meta.id}"
    """
    rgi main --input_sequence $fasta --output_file ${meta.id}_rgi --input_type contig --clean

    echo \$(rgi main --version 2>&1) > ${software}.version.txt
    """
}



