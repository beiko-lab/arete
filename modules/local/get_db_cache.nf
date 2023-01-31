// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

//process for acquiring cached databases
process GET_DB_CACHE {
    label 'process_medium'

    input:
        path(dbcache)

    output:
        path "VFDB_setA_pro.fas.gz", emit: vfdb
        path "CAZyDB.07312020.fa", emit: cazydb
        path "BacMet2_predicted_database.fasta.gz", emit: bacmet
        path "card.json", emit: card_json
        path "card.version.txt", emit: card_version
        path("""k2_standard_8gb_20201202"""), emit: minikraken
    
    script:
        """
        cp $dbcache/VFDB_setA_pro.fas.gz .
        cp $dbcache/CAZyDB.07312020.fa .
        cp $dbcache/BacMet2_predicted_database.fasta.gz .
        cp $dbcache/card.json .
        cp $dbcache/card.version.txt .
        cp -r $dbcache/k2_standard_8gb_20201202 .
        """
}

// process GET_ANNOTATION_DB_CACHE {
//     input:
//         path(dbcache)

//     output:
//         path "VFDB_setB_pro.fas.gz", emit: vfdb
//         path "CAZyDB.07312020.fa", emit: cazydb
//         path "BacMet2_predicted_database.fasta.gz", emit: bacmet
//         path "card.json", emit: card_json
//         path "card.version.txt", emit: card_version
    
//     script:
//         """
//         cp VFDB_setB_pro.fas.gz .
//         cp CAZyDB.07312020.fa .
//         cp BacMet2_predicted_database.fasta.gz .
//         cp card.json .
//         cp card.version.txt .
//         """
// }

// process GET_ASSEMBLY_DB_CACHE {
//     input:
//         path(dbcache)

//     output:
//         path("""k2_standard_8gb_20201202"""), emit: minikraken
    
//     script:
//         """
//         cp -r $dbcache/k2_standard_8gb_20201202 .
//         """
// }
