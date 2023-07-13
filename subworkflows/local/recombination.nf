include { GET_RECOMB_INPUT } from '../../modules/local/get_recomb_input'
include { VERTICALL_PAIRWISE } from '../../modules/local/verticall/pairwise.nf'
include { QUAST } from '../../modules/nf-core/quast/main'

workflow RECOMBINATION {

    take:
        assemblies
        poppunk_clusters

    main:

        ch_software_versions = Channel.empty()
        ch_reference_genome = params.reference_genome ? file(params.reference_genome) : []
        use_reference_genome = params.reference_genome ? true : false

        QUAST (
            assemblies.collect { meta, fasta -> fasta },
            ch_reference_genome,
            [],
            use_reference_genome,
            false
        )
        ch_software_versions = ch_software_versions.mix(QUAST.out.versions.ifEmpty(null))

        assemblies
            .collectFile(newLine: true) { item ->
                [ "${item[0].id}.txt", item[0].id + ',' + item[1].toString() ]
            }
            .collectFile(name: 'assembly_samplesheet.tsv')
            .set { assembly_samplesheet }

        GET_RECOMB_INPUT (
            QUAST.out.transposed_report,
            poppunk_clusters,
            assembly_samplesheet
        )

        GET_RECOMB_INPUT.out.csv
            .splitCsv(header: true)
            .map { row -> tuple(row.Cluster, file(row.samplesheet), row.Reference, row.reference_path) }
            .set { recomb_input }

        assemblies
            .combine (recomb_input.map { c,s,ref,reffile -> ref } | collect | map { [it] })
            .filter { meta, path, to_remove -> !(meta["id"] in to_remove) }
            .map { it[0, 1] }
            .collect { meta, path -> path }
            .set { unduplicated_assemblies }

        VERTICALL_PAIRWISE (
            recomb_input.map { c, s, r, rf -> tuple(c,s,rf) },
            assemblies.collect { meta, path -> path }
        )



    emit:
        verticall_result = VERTICALL_PAIRWISE.out.tsv
}
