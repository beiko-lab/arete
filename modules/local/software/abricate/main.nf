include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process UPDATE_ABRICATE_DB {
    tag "$db"
    label 'process_low'
    conda (params.enable_conda ? "bioconda::abricate=1.0.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/abricate:1.0.1--h1341992_0"
    } else {
        container "quay.io/biocontainers/abricate:1.0.1--h1341992_0"
    }
    input:
    val db

    output:
    path("""$db"""), emit: db
    
    script:
    """
    abricate-get_db --db $db --dbdir . --force
    abricate --setupdb --datadir .
    """
}

process ABRICATE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id + "/annotation/" + getSoftwareName(task.process) + "_" + db, publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::abricate=1.0.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/abricate:1.0.1--h1341992_0"
    } else {
        container "quay.io/biocontainers/abricate:1.0.1--h1341992_0"
    }

    input:
    tuple val(meta), path(fasta)
    path(db)
    val db_label

    output:
    tuple val(meta), path("${meta.id}_${db_label}.tsv"), emit: tsv
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    prefix = "${meta.id}"
    """
    abricate \\
        --datadir .\\
        --db $db_label \\
        --threads $task.cpus \\
        $fasta > ${meta.id}_${db_label}.tsv

    echo \$(abricate --version 2>&1) | sed 's/^.*abricate //' > ${software}.version.txt
    """
}




