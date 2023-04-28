#!/usr/bin/env python

import sys
import argparse
from pandas import read_csv


def parse_args(args=None):
    Description = "Add genome ID column to annotation results"
    Epilog = (
        "Example usage: python add_column.py <FILE_IN> <GENOME> <SKIP_ROWS> <FILE_OUT>"
    )

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input alignment file.")
    parser.add_argument("GENOME", help="Genome identifier.")
    parser.add_argument(
        "SKIP_ROWS", help="How many rows of text to skip in the original file."
    )
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def add_column(file_in, genome, skip_n_rows, file_out):
    ann_df = read_csv(file_in, delimiter="\t", skiprows=int(skip_n_rows))

    # Add genome_id column
    ann_df["genome_id"] = genome

    ann_df.to_csv(file_out, sep="\t", index=False)


def main(args=None):
    args = parse_args(args)
    add_column(args.FILE_IN, args.GENOME, args.SKIP_ROWS, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
