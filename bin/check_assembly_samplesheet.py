#!/usr/bin/env python

# Input check for post-assembly "Annotation" workflow

import os
import sys
from samplesheet_utils import parse_args, make_dir, print_error


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    sample,fna_file_path
    sample,fna_file_path

    """

    sample_mapping_dict = {}
    with open(file_in, "r") as fin:
        ## Check header
        MIN_COLS = 2
        HEADER = ["sample", "fna_file_path"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print(
                "ERROR: Please check samplesheet header -> {} != {}".format(
                    ",".join(header), ",".join(HEADER)
                )
            )
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(
                        MIN_COLS
                    ),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample, fna_file_path = lspl[: len(HEADER)]
            if sample:
                if sample.find(" ") != -1:
                    print_error("Sample entry contains spaces!", "Line", line)
            else:
                print_error("Sample entry has not been specified!", "Line", line)

            ## Check assembly fna file extension

            if fna_file_path:
                if fna_file_path.find(" ") != -1:
                    print_error("fna file contains spaces!", "Line", line)
                if (
                    not fna_file_path.endswith(".fna")
                    and not fna_file_path.endswith(".fna.gz")
                    and not fna_file_path.endswith(".fa")
                    and not fna_file_path.endswith(".fa.gz")
                ):
                    print_error(
                        "Assembly files must be one of .fa, .fna, .fa.gz, or .fna.gz.",
                        "Line",
                        line,
                    )
            """
            ## Auto-detect paired-end/single-end
            sample_info = []  ## [single_end, fastq_1, fastq_2]
            if sample and fastq_1 and fastq_2:  ## Paired-end short reads
                sample_info = ["0", fastq_1, fastq_2]
            elif sample and fastq_1 and not fastq_2:  ## Single-end short reads
                sample_info = ["1", fastq_1, fastq_2]
            else:
                print_error("Invalid combination of columns provided!", "Line", line)
            """
            sample_info = fna_file_path

            ## Create sample mapping dictionary = { sample: path }
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = sample_info
            # TODO Come back to this conditional; does it make sense for multiple assemblies to one sample?
            else:
                if sample_info in sample_mapping_dict[sample]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    # sample_mapping_dict[sample].append(sample_info)
                    print_error(
                        "Mutliple files associated to one sample!", "Line", line
                    )

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["sample", "path"]) + "\n")
            for sample in sorted(sample_mapping_dict.keys()):
                fout.write(
                    ",".join(
                        ["{}".format(sample), "{}".format(sample_mapping_dict[sample])]
                    )
                    + "\n"
                )
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
