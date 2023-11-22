workflow RSPR_INPUT_CHECK {
    take:
    samplesheet

    main:
    samplesheet
        .set { trees }

    emit:
    trees
}

def get_sample_info_rspr(row) {

    def array = []
    if (!file(row).exists()) {
        print("***")
        print(row)
        print("***")
        exit 1, "ERROR: Please check input samplesheet -> Tree file does not exist!\n${row.path}"
    }
    array = [ file(row)]

    return array
}
