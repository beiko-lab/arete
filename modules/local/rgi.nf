// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process UPDATE_RGI_DB {
    tag "CARD"
    label 'process_low'
    label 'error_retry_delay'

    output:
    path "card.json", emit: card_json
    path "card.version.txt", emit: card_version

    script:
    """
    curl https://card.mcmaster.ca/latest/data --output card.tar.bz2
    tar xvf card.tar.bz2
    python -c "import json;fh = open('card.json');card=json.load(fh);print(card['_version'])" > card.version.txt
    """

    stub:
    """
    touch card.json
    touch card.version.txt
    """
}

process RGI {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::rgi=6.0.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/rgi:6.0.2--pyha8f3691_0"
    } else {
        container "quay.io/biocontainers/rgi:6.0.2--pyha8f3691_0"
    }

    input:

    tuple val(meta), path(fasta)
    path card_db

    output:
    tuple val(meta), path("${meta.id}_rgi.txt"), emit: tsv
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    rgi load -i $card_db --local
    rgi main --local --input_sequence $fasta --output_file ${meta.id}_rgi --input_type contig --clean --include_loose

    echo \$(rgi main --version 2>&1) > ${software}.version.txt
    """
    stub:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    touch ${prefix}_rgi.txt
    echo \$(rgi main --version 2>&1) > ${software}.version.txt
    """
}
