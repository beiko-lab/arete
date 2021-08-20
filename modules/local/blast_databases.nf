process GET_CAZYDB {
    output:
    path "CAZyDB.07312020.fa", emit: cazydb

    script:
    """
    curl https://bcb.unl.edu/dbCAN2/download/CAZyDB.07312020.fa --output CAZyDB.07312020.fa
    """
}
