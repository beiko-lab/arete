process KRAKEN2_DB {
    label 'process_medium'
    label 'error_retry_delay'

    output:
    path("""k2_standard_8gb_20201202"""), emit: minikraken

    script:
    """
    curl https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20201202.tar.gz --output k2_standard_8gb_20201202.tar.gz
    mkdir -p k2_standard_8gb_20201202
    tar xvf k2_standard_8gb_20201202.tar.gz -C k2_standard_8gb_20201202
    """
}
