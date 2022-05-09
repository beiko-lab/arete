process GET_CAZYDB {
    output:
    path "CAZyDB.07312020.fa", emit: cazydb

    script:
    """
    curl https://bcb.unl.edu/dbCAN2/download/CAZyDB.07312020.fa --output CAZyDB.07312020.fa
    """
}

process GET_VFDB{
    output:
    path "VFDB_setB_pro.fas.gz", emit: vfdb

    script:
    """
    curl http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz --output VFDB_setB_pro.fas.gz
    """
    stub:
    """
    touch VFDB_setB_pro.fas.gz
    """
}

process GET_BACMET{
    output:
    path "BacMet2_predicted_database.fasta.gz", emit: bacmet

    script:
    """
    curl http://bacmet.biomedicine.gu.se/download/BacMet2_predicted_database.fasta.gz --output BacMet2_predicted_database.fasta.gz
    """
    stub:
    """
    touch BacMet2_predicted_database.fasta.gz
    """
}
