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
    touch VFDB_setA_pro.fas.gz
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

process GET_ICEBERG {
    label 'process_low'
    label 'error_retry_delay'

    output:
    path "ICE_aa_experimental_reformatted.fas", emit: iceberg

    script:
    """
    curl https://bioinfo-mml.sjtu.edu.cn/ICEberg2/download/ICE_aa_experimental.fas --output ICE_aa_experimental.fas
    sed -E 's/>(ICEberg\|[0-9]+)\s+(gi.*\|)\s+(.*)\s+\[(.*)\]/>\1_\3_\2_[\4]/' <ICE_aa_experimental.fas | \
    tr ' ' '_' > ICE_aa_experimental_reformatted.fas
    """
    stub:
    """
    touch ICE_aa_experimental_reformatted.fas
    """
}
