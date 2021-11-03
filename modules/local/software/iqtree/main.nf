include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process IQTREE {

    tag "roary_core"
    label 'process_high'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"phylo/iqtree") }
    
    conda (params.enable_conda ? "bioconda::iqtree=2.0.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/iqtree:2.0.3--h176a8bc_1"
    } else {
        container "quay.io/biocontainers/iqtree:2.0.3--h176a8bc_1"
    }

    input:
    path var_aln 
    path base_freq

    output:
    path "core_genome_phylo*", emit: tree
    path "*.version.txt", emit: version
    
    script:
    def software = getSoftwareName(task.process)
    """
    iqtree -fconst \$(cat $base_freq) -s $var_aln -pre core_genome_phylo -nt $task.cpus -bb 10000 -abayes -nm 1500
    echo \$(iqtree --version 2>&1) | grep "version" | sed 's/IQ-TREE multicore version //' | sed 's/ for .*//' > ${software}.version.txt
    """
}

