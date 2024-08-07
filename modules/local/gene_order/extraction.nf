process EXTRACTION {
    tag "extraction"
    label 'process_medium'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/jtllab/gene-order-workflow' :
        'docker.io/jtllab/gene-order-workflow' }"

    input:
      path input_file_path
      path extract_path
      path gbk_path
      val num_neighbors
      val percent_cutoff
      val label_cols

    output:
      path "fasta", emit: fasta_path
      path "JSON", emit: json_path
      path "neighborhood_indices.json", emit: neighborhood_indices
      path "Neighborhoods_Extraction_Summary.txt", emit: neighborhood_summary

    // This script is bundled with the pipeline, in nf-core/geneorderanalysis/bin
    script:
    """
    extraction.py \\
        -i $input_file_path \\
        -x $extract_path \\
        -g $gbk_path \\
        -o . \\
        -n $num_neighbors \\
        -p $percent_cutoff \\
        -c $label_cols
    """
}
