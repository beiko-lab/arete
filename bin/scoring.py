#!/usr/bin/env python
"""
Script for scoring utilities.
"""

import os
import re
import polars as pl
from utils import get_filename


def get_g1_g1_dict(blast_path):
    """
    Returns a dictionary of g1_g1 bitscores from BLAST output files, where keys are contig IDs (taken from the
    query_id column) and values are the corresponding bitscore value.
    """
    filepaths = os.listdir(blast_path)
    regex = r'^\w+\.\w+\.dmnd_\w+\.\w+\.txt$'
    g1_g1_dict = dict()
    for filename in filepaths:
        file_dict = dict()
        if re.search(regex, filename):
            with open(blast_path + '/' + filename, 'r') as infile:
                for line in infile.readlines():
                    tokens = line.split('\t')
                    if tokens[0] == tokens[1]:
                        g1_g1_dict[tokens[0]] = float(tokens[3].strip('\n'))

    return g1_g1_dict


def normalize_bitscores(blast_df, g1_g1_dict):
    """
    Adds normalized bitscore column to BLAST results dataframe given dataframe contains data from BLAST textfiles
    in output format 6.
    """
    blast_df = blast_df.with_column(pl.struct(['bitscore', 'query_id'])
                                    .apply(lambda row: row['bitscore'] / g1_g1_dict[row['query_id']])
                                    .alias('normalized_bitscore'))

    return blast_df


def get_normalized_bitscores(BLAST_df_dict, blast_path):
    """
    Driver method to obtain normalized bitscores for all BLAST output files in dataset.
    """
    # Keep track of contig_ids and their normalized bitscore
    g1_g1_dict = get_g1_g1_dict(blast_path)

    final_BLAST_dict = dict()
    for gene_subdir, df_dict in BLAST_df_dict.items():
        blast_files_dict = {}
        for filename, blast_df in df_dict.items():
            print("Getting normalized bitscores for gene {a} in genome {f}...".format(f=filename, a=gene_subdir))
            final_df = normalize_bitscores(blast_df, g1_g1_dict)
            blast_files_dict[filename] = final_df

        # Store data
        gene = get_filename(gene_subdir).split('.')[0]
        final_BLAST_dict[gene] = blast_files_dict

    return final_BLAST_dict
