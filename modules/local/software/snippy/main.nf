include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SNIPPY_CORE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"phylo/snippy-core") }

    conda (params.enable_conda ? "bioconda::abricate=1.0.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/abricate:1.0.1-0"
    } else {
        container "quay.io/biocontainers/abricate:1.0.1-0"
    }

    input:
    path(snippy_folders)
    reference_genome

    output:
    path("snippy_core/*.full.aln")
    path("snippy_core/*.aln")

    
    script:
    """
    snippy-core --prefix snippy_core --ref $reference_genome $snippy_folders
    echo \$(snippy --version 2>&1) | sed 's/^.*snippy //' > ${software}.version.txt
    """
}

process SNIPPY {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id + "/phylo/" + getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::snippy=4.6.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/snippy:4.6.0-1"
    } else {
        container "quay.io/biocontainers/snippy:4.6.0-1"
    }

    input:
    tuple val(meta), path(reads)
    path(reference_genome)

    output:
    path("""${meta.id}_snp_calling"""), emit: snippy_folder

    script:
    def software = getSoftwareName(task.process)
    """
    snippy --cpus $task.cpus --outdir ${meta.id}_snp_calling --ref $reference_genome --cleanup --R1 $reads[0] --R2 $reads[1] --mincov 15
    """
}

