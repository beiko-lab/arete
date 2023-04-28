#!/usr/bin/env python

import sys
import csv
import argparse


def parse_args(args=None):
    Description = "Add genome ID column to annotation results"
    Epilog = "Example usage: python add_column.py <FILE_IN> <GENOME> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input alignment file.")
    parser.add_argument("GENOME", help="Genome identifier.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def add_column(file_in, genome, file_out):
    # Add genome_id column to annotation file
    with open(file_in, "r") as r_csvfile:
        with open(file_out, "w", newline="") as w_csvfile:
            dict_reader = csv.DictReader(r_csvfile, delimiter="\t")

            # Add genome_id column
            fieldnames = dict_reader.fieldnames + ["genome_id"]
            writer_csv = csv.DictWriter(w_csvfile, fieldnames, delimiter="\t")
            writer_csv.writeheader()

            for row in dict_reader:
                row["genome_id"] = genome
                writer_csv.writerow(row)


def main(args=None):
    args = parse_args(args)
    add_column(args.FILE_IN, args.GENOME, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
