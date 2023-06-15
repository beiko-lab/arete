process CLUSTERING {
    tag "clustering"
    label 'process_high'
    label 'process_high_memory'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://jtl-lab-dev/bioinf-workflows/gene-order-workflow' :
        'jtllab/gene-order-workflow' }"

    input:
      path blast_path
      path faa_path
      path fasta_path
      val num_neighbors
      val inflation
      val epsilon
      val minpts

    output:
      path "clustering", emit: cluster_path

    script:
    """
    clustering.py $faa_path $fasta_path $blast_path . -n $num_neighbors -i $inflation -e $epsilon -m $minpts
    """
}
