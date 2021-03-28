include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MOB_INIT {
    tag "MOB_INIT"
    label 'process_medium'
    conda (params.enable_conda ? "bioconda::mob_suite=3.0.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mob_suite:3.0.1-0"
    } else {
        container "quay.io/biocontainers/mob_suite:3.0.1-0"
    }

    output:
    path "mob_suite.version.txt", emit: version

    script:
    """
    mob_init 
    mob_init --version > mob_suite.version.txt
    """
}

process MOB_RECON {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id + "/annotation/" + getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::mob_suite=3.0.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mob_suite:3.0.1-0"
    } else {
        container "quay.io/biocontainers/mob_suite:3.0.1-0"
    }

    input:
    tuple val(meta), path(fasta)
    path version

    output:
    tuple val(meta), path("""${meta.id}_mob_recon"""), emit: mob_predictions

    script:
    def software = getSoftwareName(task.process)
    """
    mob_recon --infile $fasta --num_threads $task.cpus \\
    --sample_id ${meta.id} --unicycler_contigs --run_typer \\
    --outdir ${meta.id}_mob_recon
    """
}

