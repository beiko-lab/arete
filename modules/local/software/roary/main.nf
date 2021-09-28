// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/*
If i cant figure out how to gather prokka results.
process GATHER_PROKKA{
    tag 'GATHER_PROKKA_GFF'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
    
    input:

}
*/
process ROARY {
    tag "ROARY_PANGENOME"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::roary=3.13.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/3.13.0--pl526h516909a_0"
    } else {
        container "quay.io/biocontainers/roary:3.13.0--pl526h516909a_0"
    }

    input:
    // TODO check input
    path gff

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    path("${prefix}/summary_statistics.txt"), emit: roary_summary
    path("${prefix}/gene_presence_absence.csv"), emit: roary_presence_absence_csv
    path("${prefix}/gene_presence_absence.Rtab"), emit: roary_presence_absence_rtab
    path("${prefix}/pan_genome_reference.fa"), emit: roary_pan_genome_ref
    path("${prefix}/accessory_binary_genes.fa.newick"), emit: roary_accessory_binary_genes
    path("${prefix}/accessory_graph.dot"), emit: roary_accessory_graph
    path("${prefix}/core_accessory_graph.dot"), emit: roary_core_accessory_graph
    path("${prefix}/clustered_proteins"), emit: roary_clustered_proteins
    path("${prefix}/core_gene_alignment.aln"), emit: roary_core_gene_alignment
    
    
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    prefix = "roary_pangenome_results"
    """
    roary -e -n -p $task.cpus -f $prefix $gff

    echo \$(roary -w 2>&1) | sed 's/^.*roary //' > ${software}.version.txt
    """
}
