process VERTICALL_PAIRWISE {
    tag "$cluster"
    label 'process_high'

 conda "bioconda::verticall=0.4.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/verticall:0.4.1--pyhdfd78af_0':
        'quay.io/biocontainers/verticall:0.4.1--pyhdfd78af_0' }"

    input:
    tuple val(cluster), path(assemblies), val(reference)
    path assembly_files

    output:
    tuple val(cluster), path("verticall_${prefix}.tsv"), emit: tsv
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${cluster}"
    """
    mkdir assemblies_dir

    for s in \$(cut -f2 -d"\t" $assemblies); do cp \$s assemblies_dir/; done

    verticall \\
        pairwise \\
        $args \\
        -i assemblies_dir/ \\
        --threads $task.cpus \\
        -o verticall_${prefix}.tsv \\
        -r assemblies_dir/$reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verticall: \$(echo \$(verticall --version 2>&1) | sed "s/^Verticall v//g")
    END_VERSIONS
    """
    stub:
    prefix = task.ext.prefix ?: "${cluster}"
    """
    touch verticall_${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verticall: \$(echo \$(verticall --version 2>&1) | sed "s/^Verticall v//g")
    END_VERSIONS
    """
}
