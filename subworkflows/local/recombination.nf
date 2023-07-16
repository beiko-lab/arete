include { GET_RECOMB_INPUT } from '../../modules/local/get_recomb_input'
include { VERTICALL_PAIRWISE } from '../../modules/local/verticall/pairwise.nf'
include { SKA2 } from '../../modules/local/ska'
include { QUAST } from '../../modules/nf-core/quast/main'
include { GUBBINS } from '../../modules/nf-core/gubbins/main'

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

        if (params.run_gubbins) {
            SKA2 (
                recomb_input.map { c, s, r, rf -> tuple(c,s,rf) },
                assemblies.collect { meta, path -> path }
            )

            GUBBINS (
                SKA2.out.aln
            )
        }
        if (params.run_verticall) {
            VERTICALL_PAIRWISE (
                recomb_input.map { c, s, r, rf -> tuple(c,s,rf) },
                assemblies.collect { meta, path -> path }
            )
        }

        gubbins_result = params.run_gubbins ? GUBBINS.out.stats : []
        verticall_result = params.run_verticall ? VERTICALL_PAIRWISE.out.tsv : []

    emit:
        verticall_result
        gubbins_result
}
