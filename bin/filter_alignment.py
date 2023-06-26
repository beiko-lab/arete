#!/usr/bin/env python

import os
import sys
import errno
import argparse
from pandas import read_csv, read_table


def parse_args(args=None):
    Description = "Filter alignment results."
    Epilog = "Example usage: python filter_alignment.py <FILE_IN> <GENOME> <HEADER> <PIDENT> <QCOVER> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input alignment file.")
    parser.add_argument("GENOME", help="Genome identifier.")
    parser.add_argument("HEADER", help="Header for the output file.")
    parser.add_argument(
        "PIDENT", help="Minimum match identity percentage for filtering."
    )
    parser.add_argument("QCOVER", help="Minimum coverage of each match for filtering.")
    parser.add_argument("FILE_OUT", help="Output file.")
    parser.add_argument("--parse_rgi", dest="parse_rgi", action="store_true")
    return parser.parse_args(args)


def filter_rgi(file_in, genome_id, min_pident, min_qcover, file_out):
    rgi_df = read_table(file_in)

    rgi_df.columns = rgi_df.columns.str.replace(" ", "_")
    min_qcover_scaled = float(min_qcover) * 100
    rgi_sum = rgi_df[
        (rgi_df["Best_Identities"] > int(min_pident))
        & (rgi_df["Percentage_Length_of_Reference_Sequence"] > min_qcover_scaled)
    ]
    rgi_sum = rgi_sum[["Contig", "Best_Hit_ARO", "Cut_Off", "genome_id"]].rename(
        columns={"Contig": "orf", "Best_Hit_ARO": "AMR", "Cut_Off": "rgi_cutoff"}
    )
    rgi_sum["orf"] = rgi_sum["orf"].str.rsplit("_", n=1).str.get(0)

    rgi_sum["genome_id"] = genome_id

    rgi_sum.to_csv(file_out, sep="\t", index=False)


def filter_alignment(file_in, genome_id, header, min_pident, min_qcover, file_out):
    df = read_csv(file_in, sep="\t", names=header.split(" "))

    df["genome_id"] = genome_id
    df["qcover"] = df["length"] / df["slen"]
    # Filters out every sequence that is below pident identity
    # and qcov coverage
    filtered_df = df[
        (df["pident"] >= int(min_pident)) & (df["qcover"] >= float(min_qcover))
    ]

    filtered_df = (
        filtered_df.sort_values("pident", ascending=False)
        .drop_duplicates(["qseqid", "full_qseq"])
        .sort_values(["qseqid", "pident", "mismatch"], ascending=[True, False, False])
    )

    filtered_df.to_csv(file_out, sep="\t", index=False)


def main(args=None):
    args = parse_args(args)
    if args.parse_rgi:
        filter_rgi(args.FILE_IN, args.GENOME, args.PIDENT, args.QCOVER, args.FILE_OUT)
    else:
        filter_alignment(
            args.FILE_IN,
            args.GENOME,
            args.HEADER,
            args.PIDENT,
            args.QCOVER,
            args.FILE_OUT,
        )


if __name__ == "__main__":
    sys.exit(main())
