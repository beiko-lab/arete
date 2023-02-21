process GET_CAZYDB {
    label 'process_low'
    label 'error_retry_delay'

    output:
    path "CAZyDB.07312020.fa", emit: cazydb

    script:
    """
    curl https://bcb.unl.edu/dbCAN2/download/CAZyDB.07312020.fa --output CAZyDB.07312020.fa
    """
}

process GET_VFDB{
    label 'process_low'
    label 'error_retry_delay'

    output:
    path "VFDB_setA_pro.fas.gz", emit: vfdb

    script:
    """
    curl http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz --output VFDB_setA_pro.fas.gz
    """
    stub:
    """
    touch VFDB_setB_pro.fas.gz
    """
}

process GET_BACMET{
    label 'process_low'
    label 'error_retry_delay'

    output:
    path "BacMet2_EXP_database.fasta", emit: bacmet

    script:
    """
    curl http://bacmet.biomedicine.gu.se/download/BacMet2_EXP_database.fasta --output BacMet2_EXP_database.fasta
    """
    stub:
    """
    touch BacMet2_EXP_database.fasta
    """
}
