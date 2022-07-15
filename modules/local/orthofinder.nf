process OF_PREPARE_BLAST {
    tag "$meta.id"
    label 'process_low'
    
    conda (params.enable_conda ? "bioconda::orthofinder=2.5.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orthofinder:2.5.4--hdfd78af_0':
        'quay.io/biocontainers/orthofinder:2.5.4--hdfd78af_0' }"

    input:

    tuple val(meta), path(aa_fasta)

    output:
    
    tuple val(meta), path("*.dmnd"), emit: dmnd
    tuple val(meta), path("*.fa"), emit: fa
    tuple val(meta), path('SpeciesIDs.txt'), emit: ids
    
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    orthofinder \\
        -f $aa_fasta \\
        $args \\
        -t $task.cpus \\
        -op
        -o orthofinder_results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthofinder: \$(echo \$(orthofinder -h | head -2 | tail -1 | grep -o '[0-9]+.[0-9]+.[0-9]+'))
    END_VERSIONS
    """
}

process OF_ORGANISE_BLAST{
    tag "meta.id"
    label "process_low"

    input:
    path(files_dir)

    output:
    

    tuple 
}

process ORTHOFINDER_RUN {
    tag "$meta.id"
    label 'process_high'
    conda (params.enable_conda ? "bioconda::orthofinder=2.5.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orthofinder:2.5.4--hdfd78af_0':
        'quay.io/biocontainers/orthofinder:2.5.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(blast)
    path(speciesIDs)

    output:
    path('results/*'), emit: ortho_results

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    orthofinder \\
        -b $blast \\
        $args \\
        -t $task.cpus \\
        
    """


}
