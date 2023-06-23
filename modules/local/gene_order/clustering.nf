process CLUSTERING {
    tag "clustering"
    label 'process_high'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://jtl-lab-dev/bioinf-workflows/gene-order-workflow' :
        'jtllab/gene-order-workflow' }"

    input:
      path faa_path
      path fasta_path
      path blast_path
      path json_path
      path neighborhood_indices
      val num_neighbors
      val inflation
      val epsilon
      val minpts
      path html_template
      val plot_clustering

    output:
      path "clustering", emit: cluster_path
      path "clustering_summary.csv", emit: summary
      path "JSON", emit: json

    script:
    def plot_html = plot_clustering ? "--plot_clustering" : ""
    """
    clustering.py \\
        -a $faa_path \\
        -f $fasta_path \\
        -b $blast_path \\
        -o . \\
        -n $num_neighbors \\
        -i $inflation \\
        -e $epsilon \\
        -m $minpts \\
        $plot_html
    """
}
