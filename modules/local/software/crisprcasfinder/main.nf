include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CRISPRCASFINDER {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:meta.id + "/annotation/crisprcasfinder", publish_id:meta.id) }

    container "https://crisprcas.i2bc.paris-saclay.fr/Home/DownloadFile?filename=CrisprCasFinder.simg"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple path("*.fna"), path("*.json"), emit: raw_results
    path("*.tsv"), emit: summary
    path "*.version.txt", emit: version

    script:
    """
    perl /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl -so /usr/local/CRISPRCasFinder/sel392v2.so -cf /usr/local/CRISPRCasFinder/CasFinder-2.0.3 -drpt /usr/local/CRISPRCasFinder/supplementary_files/repeatDirection.tsv -rpts /usr/local/CRISPRCasFinder/supplementary_files/Repeat_List.csv -cas -def G -gff -LOG -ccvr -gscf -meta -out ${meta.id}_out -in $fasta

    mv out/rawCRISPRs.fna ${meta.id}_rawCRISPRs.fna
    mv out/rawCas.fna ${meta.id}_rawCas.fna
    mv out/result.json ${meta.id}_result.json
    for tsv in out/TSV/*.tsv; mv \$tsv ${meta.id}_\$(basename \$tsv); done

    echo \$(perl /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl --version 2>&1) | sed 's/.* version //' | sed 's/,//' > crisprcasfinder.version.txt
    """
}

