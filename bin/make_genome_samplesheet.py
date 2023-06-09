#!/usr/bin/env python

"""
Script to generate CSV file containing the name of and path to every unique genome assembly from the provided
whole genome assembly files. Needed to apply nf-core DIAMOND MAKEDB to them within Nextflow workflow.
"""

import argparse
import itertools
import csv
import sys
import os

from utils import get_full_filepaths, remove_files, get_filename, check_output_path


def parse_args(args=None):
    Description = "Create CSV file listing each unique genome file's name and full path in a (name, path) col format."
    Epilog = "Example usage: python make_genome_samplesheet.py <ASSEMBLY_PATH> <OUTPUT_PATH>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('ASSEMBLY_PATH', metavar='asm_path', type=str,
                        help='Path to directory containing whole genome assembly files (e.g. .faa).')
    parser.add_argument('OUTPUT_PATH', metavar='output_path', type=str, help='Path to output directory where '
                                                                             'CSV file will be outputted.')
    return parser.parse_args(args)


def make_genomes_csv(assembly_path, output_path):
    """
    Driver script to obtain CSV listing all needed genome file details as needed for NF-CORE DIAMOND MAKEDB downstream.
    """
    check_output_path(output_path)

    # Get list of full file paths to genome assemblies
    full_genome_paths = [os.path.abspath(assembly_path) + '/' + file for file in os.listdir(assembly_path)]

    # Initialize column headers, row data
    fields = ['meta', 'file_path']
    row_data = []
    header = False

    # Get name and path for every whole genome assembly file
    for genome_file_path in full_genome_paths:
        row = [genome_file_path.split('/')[-1], genome_file_path]
        row_data.append(row)

    # Write data to csv file
    with open(output_path + '/genome_paths.csv', 'w') as csv_file:
        csv_writer = csv.writer(csv_file)
        if not header:
            csv_writer.writerow(fields)
            header = True
        csv_writer.writerows(row_data)


def main(args=None):
    args = parse_args(args)
    make_genomes_csv(args.ASSEMBLY_PATH, args.OUTPUT_PATH)


if __name__ == '__main__':
    sys.exit(main())
