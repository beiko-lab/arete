//
// MODULE: Installed directly from nf-core/modules
//
include { PROKKA                } from '../../modules/nf-core/prokka/main'
include { BAKTA_BAKTA as BAKTA } from '../../modules/nf-core/bakta/bakta/main'
include { BAKTA_BAKTADBDOWNLOAD as BAKTADBDOWNLOAD } from '../../modules/nf-core/bakta/baktadbdownload/main'
include { GET_CAZYDB;
          GET_VFDB;
          GET_BACMET;
          GET_ICEBERG } from '../../modules/local/blast_databases.nf'
include { ADD_GENOME_COLUMN as PROKKA_ADD_COLUMN;
          ADD_GENOME_COLUMN as BAKTA_ADD_COLUMN;
          ADD_GENOME_COLUMN as RGI_ADD_COLUMN } from '../../modules/local/add_genome_column'
include { DIAMOND_MAKEDB as DIAMOND_MAKE_CAZY;
          DIAMOND_MAKEDB as DIAMOND_MAKE_VFDB;
          DIAMOND_MAKEDB as DIAMOND_MAKE_BACMET;
          DIAMOND_MAKEDB as DIAMOND_MAKE_ICEBERG } from '../../modules/nf-core/diamond/makedb/main'
include { DIAMOND_BLASTX as DIAMOND_BLAST_CAZY;
          DIAMOND_BLASTX as DIAMOND_BLAST_VFDB;
          DIAMOND_BLASTX as DIAMOND_BLAST_BACMET;
          DIAMOND_BLASTX as DIAMOND_BLAST_ICEBERG } from '../../modules/nf-core/diamond/blastx/main'
//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../../modules/local/get_software_versions'
include { RGI;
          UPDATE_RGI_DB } from '../../modules/local/rgi'
include { MOB_RECON } from '../../modules/local/mobsuite'
include { ISLANDPATH } from '../../modules/local/islandpath/main'
include { VIBRANT_DOWNLOADDB } from '../../modules/local/vibrant/downloadb/main.nf'
include { VIBRANT_VIBRANTRUN } from '../../modules/local/vibrant/vibrantrun/main.nf'
include { INTEGRON_FINDER } from '../../modules/local/integronfinder/main.nf'
include { CONCAT_OUTPUT as CONCAT_PROKKA;
          CONCAT_OUTPUT as CONCAT_BAKTA;
          CONCAT_OUTPUT as CONCAT_RGI;
          CONCAT_OUTPUT as CONCAT_MOBSUITE;
          CONCAT_OUTPUT as CONCAT_ISLANDS;
          CONCAT_OUTPUT as CONCAT_INTEGRONS } from '../../modules/local/concat_output.nf'

//
// SUBWORKFLOWS
//
include { FILTER_ALIGNMENT as CAZY_FILTER;
          FILTER_ALIGNMENT as VFDB_FILTER;
          FILTER_ALIGNMENT as BACMET_FILTER;
          FILTER_ALIGNMENT as ICEBERG_FILTER } from './concatenate_matches'



workflow ANNOTATE_ASSEMBLIES {
    take:
        assemblies
        bakta_db
        vfdb_cache
        cazydb_cache
        bacmet_cache
        icebergdb_cache
        card_json_cache
        card_version_cache


    main:

        ch_multiqc_files = Channel.empty()
        ch_software_versions = Channel.empty()

        min_pident = params.min_pident
        min_qcover = params.min_qcover
        /*
        * SUBWORKFLOW: Read in samplesheet, validate and stage input files
        */

        /*
        * Load in the databases. Check if they were cached, otherwise run the processes that get them
        */

        // Note: I hate how inelegant this block is. Works for now, but consider looking for a more elegant nextflow pattern
        /*
        * Load BLAST databases
        */
        if (vfdb_cache){
            vfdb_cache.set { ch_vfdb }
        }
        else{
            GET_VFDB()
            GET_VFDB.out.vfdb.set { ch_vfdb }
        }
        if (!params.light) {
            if(bacmet_cache){
                bacmet_cache.set { ch_bacmet_db }
            }
            else{
                GET_BACMET()
                GET_BACMET.out.bacmet.set { ch_bacmet_db }
            }
            if (cazydb_cache){
                cazydb_cache.set{ ch_cazy_db }
            }
            else{
                GET_CAZYDB()
                GET_CAZYDB.out.cazydb.set { ch_cazy_db }
            }
            if (icebergdb_cache){
                icebergdb_cache.set{ ch_iceberg_db }
            }
            else{
                GET_ICEBERG()
                GET_ICEBERG.out.iceberg.set { ch_iceberg_db }
            }
        }
        /*
        * Load RGI for AMR annotation
        */
        if (card_json_cache){
            card_json_cache.set { ch_card_json }
            ch_software_versions = ch_software_versions.mix(card_version_cache)
        }
        else{
            UPDATE_RGI_DB()
            UPDATE_RGI_DB.out.card_json.set { ch_card_json }
            ch_software_versions = ch_software_versions.mix(UPDATE_RGI_DB.out.card_version.ifEmpty(null))
        }

        /*
        * Run gene finding software (Prokka or Bakta)
        */
        if (params.use_prokka) {

            PROKKA (
            assemblies,
            [],
            []
            ) //Assembly, protein file, pre-trained prodigal
            ch_software_versions = ch_software_versions.mix(PROKKA.out.versions.first().ifEmpty(null))
            ch_ffn_files = PROKKA.out.ffn
            ch_gff_files = PROKKA.out.gff
            ch_gbk_files = PROKKA.out.gbk
            ch_tsv_files = PROKKA.out.tsv
            ch_multiqc_files = ch_multiqc_files.mix(PROKKA.out.txt.collect{it[1]}.ifEmpty([]))

            PROKKA_ADD_COLUMN(
                ch_tsv_files,
                "PROKKA",
                0
            )

            PROKKA_ADD_COLUMN.out.txt
                .collect{ id, path -> path }
                .set{ prokka_tsvs }

            CONCAT_PROKKA(prokka_tsvs, "PROKKA", 1)
        }
        else {

            if (bakta_db){
                BAKTA(assemblies, bakta_db, [], [])
            } else {
                BAKTADBDOWNLOAD()
                BAKTADBDOWNLOAD.out.db.set { bakta_db }
                BAKTA(assemblies, bakta_db, [], [])
            }

            ch_software_versions = ch_software_versions.mix(BAKTA.out.versions.first().ifEmpty(null))
            ch_ffn_files = BAKTA.out.ffn
            ch_gff_files = BAKTA.out.gff
            ch_gbk_files = BAKTA.out.gbff
            ch_tsv_files = BAKTA.out.tsv

            BAKTA_ADD_COLUMN(
                ch_tsv_files,
                "BAKTA",
                2
            )

            BAKTA_ADD_COLUMN.out.txt
                .collect{ id, path -> path }
                .set{ bakta_tsvs }

            CONCAT_BAKTA(bakta_tsvs, "BAKTA", 1)
        }

        if (!params.skip_vibrant){
            if (params.vibrant_db){
                VIBRANT_VIBRANTRUN(assemblies, file(params.vibrant_db))
            } else {
                VIBRANT_DOWNLOADDB()
                VIBRANT_DOWNLOADDB.out.db.set { vibrant_db }
                VIBRANT_VIBRANTRUN(assemblies, vibrant_db)
            }
        }

        /*
        * Run RGI
        */
        RGI(ch_ffn_files, ch_card_json)
        ch_software_versions = ch_software_versions.mix(RGI.out.version.first().ifEmpty(null))

        RGI_ADD_COLUMN(
          RGI.out.tsv,
          "RGI",
          0
        )

        RGI_ADD_COLUMN.out.txt
            .collect{ id, paths -> paths }
            .set { rgi_tsvs }

        CONCAT_RGI(rgi_tsvs, "RGI", 1)

        /*
        * Module: Mob-Suite. Database is included in singularity container
        */
        MOB_RECON(assemblies)
        ch_software_versions = ch_software_versions.mix(MOB_RECON.out.version.first().ifEmpty(null))

        MOB_RECON.out.contig_report
            .collect{ id, paths -> paths }
            .set { mobrecon_tsvs }

        CONCAT_MOBSUITE(mobrecon_tsvs, "MOBSUITE", 1)

        if (params.run_integronfinder){
            INTEGRON_FINDER(assemblies)
            ch_software_versions = ch_software_versions.mix(INTEGRON_FINDER.out.versions.first())

            INTEGRON_FINDER.out.summary
                .collect{ id, paths -> paths }
                .set{ integron_summaries }

            CONCAT_INTEGRONS(integron_summaries, "INTEGRONFINDER", 2)
        }

        ISLANDPATH(ch_gbk_files)
        ch_software_versions = ch_software_versions.mix(ISLANDPATH.out.versions.first())

        ISLANDPATH.out.gff
            .collect{ id, paths -> paths }
            .set { islandpath_gffs }

        CONCAT_ISLANDS(islandpath_gffs, "ISLANDPATH", 1)

        /*
        * Run DIAMOND blast annotation with databases
        */
        def blast_columns = "qseqid sseqid pident slen qlen length mismatch gapopen qstart qend sstart send evalue bitscore full_qseq"


        DIAMOND_MAKE_VFDB(ch_vfdb)
        DIAMOND_BLAST_VFDB(ch_ffn_files, DIAMOND_MAKE_VFDB.out.db, "txt", blast_columns)
        VFDB_FILTER(
            DIAMOND_BLAST_VFDB.out.txt,
            "VFDB",
            blast_columns,
            min_pident,
            min_qcover
        )

        if (!params.light) {
            DIAMOND_MAKE_BACMET(ch_bacmet_db)
            DIAMOND_BLAST_BACMET(ch_ffn_files, DIAMOND_MAKE_BACMET.out.db, "txt", blast_columns)
            BACMET_FILTER(
                DIAMOND_BLAST_BACMET.out.txt,
                "BACMET",
                blast_columns,
                min_pident,
                min_qcover
            )

            DIAMOND_MAKE_CAZY(ch_cazy_db)
            DIAMOND_BLAST_CAZY(ch_ffn_files, DIAMOND_MAKE_CAZY.out.db, "txt", blast_columns)
            CAZY_FILTER(
                DIAMOND_BLAST_CAZY.out.txt,
                "CAZY",
                blast_columns,
                min_pident,
                min_qcover
            )

            DIAMOND_MAKE_ICEBERG(ch_iceberg_db)
            DIAMOND_BLAST_ICEBERG(ch_ffn_files, DIAMOND_MAKE_ICEBERG.out.db, "txt", blast_columns)
            ICEBERG_FILTER(
                DIAMOND_BLAST_ICEBERG.out.txt,
                "ICEBERG",
                blast_columns,
                min_pident,
                min_qcover
            )
        }

        ch_software_versions = ch_software_versions.mix(DIAMOND_MAKE_VFDB.out.versions.ifEmpty(null))

    emit:
        annotation_software = ch_software_versions
        multiqc = ch_multiqc_files
        gff = ch_gff_files

}
