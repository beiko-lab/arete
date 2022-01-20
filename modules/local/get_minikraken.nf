// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process KRAKEN2_DB {
    tag "minikraken"
    label 'process_high'

    output:
    path("""k2_standard_8gb_20201202"""), emit: minikraken

    script:
    """
    curl https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20201202.tar.gz --output k2_standard_8gb_20201202.tar.gz
    mkdir -p k2_standard_8gb_20201202
    tar xvf k2_standard_8gb_20201202.tar.gz -C k2_standard_8gb_20201202
    """
}
