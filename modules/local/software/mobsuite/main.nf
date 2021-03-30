include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MOB_RECON {
    tag "$meta.id"
    label 'process_low'
    stageInMode 'copy'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id + "/annotation/mob_recon", publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::mob_suite=3.0.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mob_suite:3.0.1-0"
    } else {
        container "finlaymaguire/mob-suite:latest"
        //container "quay.io/biocontainers/mob_suite:3.0.1--py_0"
    }

    input:
    tuple val(meta), path(fasta)
    path mob_db

    output:
    tuple val(meta), path("""${meta.id}_mob_recon"""), emit: mob_predictions
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    mob_recon --infile $fasta --num_threads $task.cpus \\
    --sample_id ${meta.id} --unicycler_contigs --run_typer \\
    --outdir ${meta.id}_mob_recon --database_directory $mob_db  \\
    --run_overhang
    mob_recon --version > ${software}.version.txt
    """
}

process MOB_INIT {
    tag "MOB_INIT"
    label 'process_medium'
    conda (params.enable_conda ? "bioconda::mob_suite=3.0.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mob_suite:3.0.1-0"
    } else {
        container "finlaymaguire/mob-suite:latest"
        //container "quay.io/biocontainers/mob_suite:3.0.1--py_0"
    }

    output:
    path """mob_db""", emit: mob_db
    path "*.version.txt", emit: version

    script:
    """
    mkdir -p mob_db
    mob_init -d mob_db
    mob_init --version > mob_suite.version.txt
    """
}

