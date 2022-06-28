def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

def modules = params.modules.clone()

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { PROKKA                } from '../../modules/nf-core/modules/prokka/main' addParams( options: [:] )
include { GET_CAZYDB;
          GET_VFDB;
          GET_BACMET} from '../../modules/local/blast_databases.nf'
include { DIAMOND_MAKEDB as DIAMOND_MAKE_CAZY;
          DIAMOND_MAKEDB as DIAMOND_MAKE_VFDB;
          DIAMOND_MAKEDB as DIAMOND_MAKE_BACMET } from '../../modules/nf-core/modules/diamond/makedb/main'  addParams( options: [:] )
include { DIAMOND_BLASTX as DIAMOND_BLAST_CAZY;
          DIAMOND_BLASTX as DIAMOND_BLAST_VFDB;
          DIAMOND_BLASTX as DIAMOND_BLAST_BACMET } from '../../modules/nf-core/modules/diamond/blastx/main'  addParams( options: [args:'--evalue 1e-06 --outfmt 6 qseqid sseqid pident slen qlen length mismatch gapopen qstart qend sstart send evalue bitscore full_qseq --max-target-seqs 25 --more-sensitive'] )
//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include { RGI;
          UPDATE_RGI_DB } from '../../modules/local/rgi'  addParams( options: [:] )
include { MOB_RECON } from '../../modules/local/mobsuite'  addParams( options: [:] )


// Usage pattern from nf-core/rnaseq: Empty dummy file for optional inputs
ch_dummy_input = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)


workflow ANNOTATE_ASSEMBLIES {
    take:
        assemblies
    main:
    
        //if (params.input_sample_table){ ch_input = file(params.input_sample_table) } else { exit 1, 'Input samplesheet not specified!' }

        ch_software_versions = Channel.empty()
        /*
        * SUBWORKFLOW: Read in samplesheet, validate and stage input files
        */
        //ANNOTATION_INPUT_CHECK(ch_input)

        /*
        * Module: Annotate AMR
        */
        UPDATE_RGI_DB()
        ch_software_versions = ch_software_versions.mix(UPDATE_RGI_DB.out.card_version.ifEmpty(null))
        RGI(assemblies, UPDATE_RGI_DB.out.card_json)
        ch_software_versions = ch_software_versions.mix(RGI.out.version.first().ifEmpty(null))

        /*
        * Module: Prokka
        */
        //TODO prokka is in both annotation and assembly right now...
        PROKKA (
        assemblies,
        [],
        []
        ) //Assembly, protein file, pre-trained prodigal
        ch_software_versions = ch_software_versions.mix(PROKKA.out.versions.first().ifEmpty(null))


        /*
        * Module: Mob-Suite
        */
        MOB_RECON(assemblies)
        ch_software_versions = ch_software_versions.mix(MOB_RECON.out.version.first().ifEmpty(null))

        /*
        * Module: BLAST vs CAZY, VFDB, Bacmet
        */
        GET_CAZYDB()
        GET_BACMET()
        GET_VFDB()
        DIAMOND_MAKE_CAZY(GET_CAZYDB.out.cazydb)
        ch_software_versions = ch_software_versions.mix(DIAMOND_MAKE_CAZY.out.versions.ifEmpty(null))
        DIAMOND_BLAST_CAZY(PROKKA.out.ffn, DIAMOND_MAKE_CAZY.out.db, "CAZYDB")

        DIAMOND_MAKE_VFDB(GET_VFDB.out.vfdb)
        DIAMOND_BLAST_VFDB(PROKKA.out.ffn, DIAMOND_MAKE_VFDB.out.db, "VFDB")

        DIAMOND_MAKE_BACMET(GET_BACMET.out.bacmet)
        DIAMOND_BLAST_BACMET(PROKKA.out.ffn, DIAMOND_MAKE_BACMET.out.db, "BACMET")


        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(PROKKA.out.txt.collect{it[1]}.ifEmpty([]))

    emit:
        annotation_software = ch_software_versions
        multiqc = ch_multiqc_files
        gff = PROKKA.out.gff

}