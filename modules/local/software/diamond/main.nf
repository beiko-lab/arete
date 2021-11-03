include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DIAMOND_MAKEDB {
    tag "$fasta"
    label 'process_low'
    conda (params.enable_conda ? "bioconda::diamond=2.0.8" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/diamond:2.0.8--h56fc30b_0"
    } else {
        container "quay.io/biocontainers/diamond:2.0.8--h56fc30b_0"
    }
    input:
    path fasta

    output:
    path("*.dmnd"), emit: db
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    diamond makedb --in $fasta -d ${fasta.baseName}
    echo \$(diamond --version 2>&1) | sed 's/^.*diamond version //' > ${software}.version.txt
    """
}

process DIAMOND_BLASTX {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id + "/annotation/" + getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::diamond=2.0.8" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/diamond:2.0.8--h56fc30b_0"
    } else {
        container "quay.io/biocontainers/diamond:2.0.8--h56fc30b_0"
    }

    input:
    tuple val(meta), path(orfs)
    path db
    val label

    output:
    path("${meta.id}_${label}.out6"), emit: out6

    script:
    def software = getSoftwareName(task.process)
    """
    diamond blastx --query $orfs --db $db --evalue 1e-06 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore full_qseq --max-target-seqs 25 --out ${meta.id}_${label}.out6 --threads $task.cpus --more-sensitive
    """
}



