// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MOB_RECON {
    tag "$meta.id"
    label 'process_low'
    //stageInMode 'copy'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id + "/annotation/mob_recon", publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::mob_suite=3.0.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mob_suite:3.0.3--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/mob_suite:3.0.3--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path(fasta)
    //path mob_db

    output:
    tuple val(meta), path("""${meta.id}_mob_recon"""), emit: mob_predictions
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    mob_recon --infile $fasta --num_threads $task.cpus \\
    --sample_id ${meta.id} --unicycler_contigs --run_typer \\
    --outdir ${meta.id}_mob_recon \\
    --run_overhang
    mob_recon --version > ${software}.version.txt
    """
    //--database_directory $mob_db  \\
    stub:
    def software = getSoftwareName(task.process)
    """
    mkdir ${meta.id}_mob_recon
    mob_recon --version > ${software}.version.txt
    """
}

// MOB_INIT also exists but seems to be unneeded.
