#!/usr/bin/env python

"""
Given a textfile where each row lists a gene to extract neighborhoods for, generates extraction input files from
Genbank records by fetching respective occurrences of the gene in every genome in the provided dataset.
"""

import os
import sys
import glob
import argparse
import pandas as pd
from Bio import SeqIO
from utils import get_filename


def parse_args(args=None):
    Description = "Generates required extraction data files for workflow from Genbank (GBK) files given a list of ."
    Epilog = "Example usage: python make_extraction_files.py <GBK_PATH> <TEXTFILE_PATH>"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('GBK_PATH', metavar='gbk', type=str, help='Path to directory containing GBK files. \
                                                                       Must have names corresponding to RGI files.')
    parser.add_argument('TEXTFILE_PATH', metavar='txt', type=str, help='Path to textfile listing one gene per row to \
                                                                        extract neighborhoods for from GBK files.')
    parser.add_argument('OUTPUT_PATH', metavar='txt', type=str, help='Output path to generate extraction textfiles in.')
    return parser.parse_args(args)


def get_genes_list(textfile_path):
    """
    Processes user inputted textfile into a list of unique genes to be extracted.
    """
    genes_list = []
    with open(textfile_path, 'r') as infile:
        try:
            text_file_data = infile.readlines()
        except Exception as e:
            sys.stderr.write('Error occurred attempting to process input textfile: {e}'.format(e=e))

    for line in text_file_data:
        gene = line.strip('\n')
        genes_list.append(gene)

    return genes_list


def get_file_paths(gbk_path):
    """
    Loads all GBK filepaths from user provided directory path and returns them in a list.
    Assumes that GBK files have the <.gbk> file extension.
    """
    try:
        gbk_filepaths = glob.glob(os.path.join(gbk_path, "*.gbk"))
    except FileNotFoundError:
        print("Error: there are no GBK files found in the specified directory. Please double-check the provided path.")
        sys.exit(1)

    return gbk_filepaths


def extract_records(gbk_file_path, genes_list):
    """
    Given gene to extract, searches through Genbank data for any occurrences in each respective genome.
    Outputs a dictionary where keys correspond to genome IDs and values are tab-delimited rows of data
    needed for the final extraction files.
    """
    records_list = []
    for gene in genes_list:
        for index, record in enumerate(SeqIO.parse(gbk_file_path, "genbank")):
            for feature in record.features:
                if feature.type == 'CDS' and 'pseudo' not in feature.qualifiers and 'gene' in feature.qualifiers:
                    for key, val in feature.qualifiers.items():
                        if gene in val:
                            data = [record.id, feature.location.start, feature.location.end, feature.qualifiers['gene']]
                            records_list.append(data)

    return records_list


def make_genome_extraction_data_dict(gbk_file_paths, genes_list):
    """
    For every GBK file, creates a dictionary entry consisting of a list of all gene occurrences and their associated
    contig name, start index, stop index, and the gene name.
    """
    genome_extraction_data = {}
    for gbk_file_path in gbk_file_paths:
        extraction_data = extract_records(gbk_file_path)
        genome_id = get_filename(gbk_file_path)
        genome_extraction_data[genome_id] = extraction_data

    return genome_extraction_data


def make_genome_extraction_input_file(genome_id, genome_data, output_path):
    """
    From data extracted from Genbank files for user selected genes, generates extraction input files for workflow
    with required data in a tab-delimited format.
    Example:
    '''
    Contig_Name   Start     Stop     Name
    04            211727    214876   acrB_1
    '''
    """
    with open(output_path + '/' + genome_id + '.txt', 'w') as outfile:
        outfile.write('Contig_Name\tStart\tStop\tName')
        for gene_instance in genome_data:
            for i in range(len(gene_instance)):
                outfile.write(gene_instance[i])
                if i != len(gene_instance) - 1:
                    outfile.write('\t')
                else:
                    outfile.write('\n')


def make_input_files(gbk_path, textfile_path, output_path):
    """
    Driver script for generating extraction input files for each genome from user provided textfile indicating genes
    to extract neighborhoods for.
    """
    # Process textfile to get gene names as list
    genes_to_extract = get_genes_list(textfile_path)

    # Get GBK filepaths
    gbk_file_paths = get_file_paths(gbk_path)

    # Get gene data for all genomes to write to final files
    extraction_dict = make_genome_extraction_data_dict(gbk_file_paths, genes_to_extract)

    # Generate required input files for workflow
    for genome_id, genome_extraction_data in extraction_dict.items():
        make_genome_extraction_input_file(genome_id, genome_extraction_data, output_path)


def main(args=None):
    args = parse_args(args)
    make_input_files(args.GBK_PATH, args.TEXTFILE_PATH, args.OUTPUT_PATH)


if __name__ == '__main__':
    sys.exit(main())
