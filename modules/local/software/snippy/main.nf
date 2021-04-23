include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SNIPPY_CORE {
    tag "snippy_core"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"phylo/snippy-core") }
    
    conda (params.enable_conda ? "bioconda::snippy=4.6.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/snippy:4.6.0--0"
    } else {
        container "quay.io/biocontainers/snippy:4.6.0--0"
    }

    input:
    path snippy_folders
    path reference_genome

    output:
    path "snippy_core.full.aln", emit: full_aln
    path "snippy_core.aln", emit: var_aln
    path "base_frequencies.txt", emit: base_freq
    path "*.version.txt", emit: version
    
    script:
    def software = getSoftwareName(task.process)
    """
    snippy-core --prefix snippy_core --ref $reference_genome $snippy_folders
    echo \$(snippy --version 2>&1) | sed 's/^.*snippy //' > ${software}.version.txt

    snp-sites -C snippy_core.full.aln > base_frequencies.txt
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
        container "https://depot.galaxyproject.org/singularity/snippy:4.6.0--0"
    } else {
        container "quay.io/biocontainers/snippy:4.6.0--0"
    }

    input:
    tuple val(meta), path(reads)
    path(reference_genome)

    output:
    path("""${meta.id}_snp_calling"""), emit: snippy_folder

    script:
    """
    snippy --cpus $task.cpus --outdir ${meta.id}_snp_calling --ref $reference_genome --cleanup --R1 $reads[0] --R2 $reads[1] --mincov 15
    rm ${meta.id}_snp_calling/ref.fa.fai
    """
}


process SNIPPY_CTG {
    tag "SNIPPY_OUTGROUP"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"phylo/outgroup_" + getSoftwareName(task.process)) }

    conda (params.enable_conda ? "bioconda::snippy=4.6.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/snippy:4.6.0--0"
    } else {
        container "quay.io/biocontainers/snippy:4.6.0--0"
    }

    input:
    path(outgroup_fasta)
    path(reference_genome)

    output:
    path("""outgroup_snp_calling"""), emit: snippy_folder

    script:
    """
    snippy --cpus $task.cpus --outdir outgroup_snp_calling --ref $reference_genome --cleanup --ctgs $outgroup_fasta --mincov 15
    rm outgroup_snp_calling/ref.fa.fai
    """
}



