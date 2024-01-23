process MOB_RECON {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::mob_suite=3.0.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mob_suite:3.0.3--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/mob_suite:3.0.3--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_mob_recon")                             , emit: mob_predictions
    tuple val(meta), path("${meta.id}_mob_recon/${meta.id}_contig_report.txt"), emit: contig_report
    path "versions.yml"                                                       , emit: version

    script:
    """
    mob_recon --infile $fasta --num_threads $task.cpus \\
    --sample_id ${meta.id} --unicycler_contigs --run_typer \\
    --outdir ${meta.id}_mob_recon \\
    --run_overhang

    mv ${meta.id}_mob_recon/contig_report.txt ${meta.id}_mob_recon/${meta.id}_contig_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mobsuite: \$(echo \$(mob_recon --version 2>&1) | sed 's/^.*mob_recon //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir ${meta.id}_mob_recon
    touch ${meta.id}_mob_recon/${meta.id}_contig_report.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mobsuite: \$(echo \$(mob_recon --version 2>&1) | sed 's/^.*mob_recon //; s/ .*\$//')
    END_VERSIONS
    """
}

// MOB_INIT also exists but seems to be unneeded.
