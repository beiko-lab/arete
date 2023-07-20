//process for acquiring cached databases
process GET_DB_CACHE {
    label 'process_single'

    input:
        path(dbcache)

    output:
        path "VFDB_setA_pro.fas.gz", emit: vfdb
        path "CAZyDB.07312020.fa", emit: cazydb
        path "BacMet2_EXP_database.fasta", emit: bacmet
        path "ICE_aa_experimental_reformatted.fas", emit: iceberg
        path "card.json", emit: card_json
        path "card.version.txt", emit: card_version
        path "k2_standard_8gb_20201202", emit: minikraken

    script:
        """
        cp $dbcache/VFDB_setA_pro.fas.gz .
        cp $dbcache/CAZyDB.07312020.fa .
        cp $dbcache/BacMet2_EXP_database.fasta .
        cp $dbcache/ICE_aa_experimental_reformatted.fas .
        cp $dbcache/card.json .
        cp $dbcache/card.version.txt .
        cp -r $dbcache/k2_standard_8gb_20201202 .
        """
}
