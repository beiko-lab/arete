process VIBRANT_VIBRANTRUN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::vibrant=1.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vibrant:1.2.1--hdfd78af_2':
        'quay.io/biocontainers/vibrant:1.2.1--hdfd78af_2' }"

    input:

    tuple val(meta), path(fasta)
    path db

    output:

    path "${prefix}"            , emit: vibrant_results
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    VIBRANT_run.py \\
        -i $fasta \\
        -d $PWD/$db/databases/ \\
        -m $PWD/$db/files/ \\
        -folder $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        VIBRANT: \$(echo \$(VIBRANT_run.py --version 2>&1) | sed 's/^.*VIBRANT v//;')
    END_VERSIONS
    """
    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir $prefix
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        VIBRANT: \$(echo \$(VIBRANT_run.py --version 2>&1) | sed 's/^.*VIBRANT v//;')
    END_VERSIONS
    """
}
