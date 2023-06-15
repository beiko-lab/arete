#!/usr/bin/env python

"""
Script to generate CSV file containing each unique filepair from among provided whole genome assembly files.
Needed in order to do All-vs-All BLAST of unique pairs using nf-core DIAMOND BLASTP within Nextflow workflow.
"""

import argparse
import itertools
import csv
import sys
import os

from utils import get_full_filepaths, remove_files, get_filename, check_output_path


def parse_args(args=None):
    Description = "Create CSV file listing each unique combination of genomes to BLAST in (genome1, genome2) format."
    Epilog = "Example usage: python make_genome_filepairs.py <ASSEMBLY_PATH> <OUTPUT_PATH>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('ASSEMBLY_PATH', metavar='asm_path', type=str,
                        help='Path to directory containing whole genome assembly files to BLAST (e.g. .faa).')
    parser.add_argument('OUTPUT_PATH', metavar='output_path', type=str, help='Path to output directory where '
                                                                             'CSV file will be outputted.')
    return parser.parse_args(args)


def make_filepairs_csv(assembly_path, output_path):
    """
    Driver script to obtain CSV listing all unique genome file pairs as needed for All-vs-All BLAST downstream
    using NF-CORE DIAMOND BLASTP.
    """
    check_output_path(output_path)

    # Get list of full file paths to genome assemblies
    full_genome_paths = [os.path.abspath(assembly_path) + '/' + file for file in os.listdir(assembly_path)]

    # Initialize column headers, row data
    fields = ['meta', 'fasta', 'db']
    row_data = []
    header = False

    # Get every identical combination of genomes
    for genome in full_genome_paths:
        g1 = get_filename(genome)
        db_filename = os.path.basename(genome)
        row = [g1 + '_' + g1, genome, os.path.abspath(output_path) + '/diamond/' + db_filename + '.dmnd']
        row_data.append(row)

    # Get every unique combination of genomes to BLAST against each other
    for genome_1, genome_2 in itertools.combinations(full_genome_paths, 2):
        g1 = get_filename(genome_1)
        g2 = get_filename(genome_2)
        db_filename = os.path.basename(genome_2)
        row = [g1 + '_' + g2, genome_1, os.path.abspath(output_path) + '/diamond/' + db_filename + '.dmnd']
        row_data.append(row)

    # Write data to csv file
    with open(output_path + '/genome_filepairs.csv', 'w+') as csv_file:
        csv_writer = csv.writer(csv_file)
        if not header:
            csv_writer.writerow(fields)
            header = True
        csv_writer.writerows(row_data)


def main(args=None):
    args = parse_args(args)
    make_filepairs_csv(args.ASSEMBLY_PATH, args.OUTPUT_PATH)


if __name__ == '__main__':
    sys.exit(main())
