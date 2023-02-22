#!/usr/bin/env python

import os
import sys
import errno
import argparse
from pandas import read_csv


def parse_args(args=None):
    Description = "Filter alignment results."
    Epilog = "Example usage: python filter_alignment.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input alignment file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def filter_alignment(file_in, file_out):
    df = read_csv(file_in, sep="\t")

    df["qcover"] = df["length"] / df["slen"]
    filtered_df = (
        df.sort_values("pident", ascending=False)
        .drop_duplicates(["qseqid", "full_qseq"])
        .sort_values(["qseqid", "pident", "mismatch"], ascending=[True, False, False])
    )

    filtered_df.to_csv(file_out, sep="\t", index=False)


def main(args=None):
    args = parse_args(args)
    filter_alignment(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
