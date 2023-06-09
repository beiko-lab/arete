#!/usr/bin/env python
"""
Functions to cluster neighborhoods identified using extraction.py.
"""

import argparse
import itertools
import os
import sys

import numpy as np
import pandas as pd
import polars as pl
from scipy import sparse
from scipy.cluster import hierarchy
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler

import markov_clustering as mcl

from scoring import get_normalized_bitscores
from utils import get_filename, check_output_path, get_full_filepaths, remove_files, generate_alphanumeric_string
from json_utils import order_JSON_clusters_UPGMA, make_representative_UPGMA_cluster_JSON, update_JSON_links_PI
from visualization import plot_similarity_histogram, plot_distance_histogram, plotly_pcoa, \
    plotly_dendrogram, plotly_mcl_network


def parse_args(args=None):
    Description = "Cluster extracted  gene neighborhoods to compare conservation characteristics across genomes."
    Epilog = "Example usage: python clustering.py <ASSEMBLY_PATH> <FASTA_PATH> <BLAST_PATH> <OUTPUT_PATH> \
              -n <NEIGHBORHOOD_SIZE>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('ASSEMBLY_PATH', metavar='asm_path', type=str,
                        help='Path to FA assembly files.')
    parser.add_argument('FASTA_PATH', metavar='fasta_path', type=str,
                        help='Path to directory containing neighborhood FASTA files.')
    parser.add_argument('BLAST_PATH', metavar='blast_path', type=str,
                        help='Path to directory containing BLAST text files of BLASTed neighborhoods.')
    parser.add_argument('OUTPUT_PATH', metavar='output_path', type=str, help='Path to output directory where '
                                                                             'extracted neighborhood FASTA files will'
                                                                             ' be saved.')
    parser.add_argument('-n', metavar='-neighborhood_size', type=int, default=10, help='Neighborhood window size, i.e. '
                                                                                       'number of genes considered upstream and '
                                                                                       'downstream of focal gene.')
    parser.add_argument('-i', metavar='-mcl_inf', type=int, default=2, help='Inflation hyperparameter for Markov \
                                                                                                           clustering.')
    parser.add_argument('-e', metavar='-dbscan_eps', type=float, default=0.5, help='Inflation hyperparameter for \
                                                                                                    DBSCAN clustering.')
    parser.add_argument('-m', metavar='-dbscan_min', type=int, default=5, help='Minimum samples hyperparameter for \
                                                                                                    DBSCAN clustering.')

    return parser.parse_args(args)


def make_fasta_contig_dict(fasta_path, gene):
    """
    Given FASTA neighborhood file, creates a dictionary where keys correspond to indices from 0 to N,
    and values are contig ids.
    """
    fasta_dict = dict()
    fasta_files = [fasta_file for fasta_file in os.listdir(fasta_path + '/' + gene)
                   if fasta_file.endswith('.fasta')]
    for fasta_file in fasta_files:
        with open(fasta_path + '/' + gene + '/' + fasta_file, 'r') as infile:
            data = infile.readlines()
            fasta_file_dict = dict()
            index = 0
            contig_lines = [line for line in data if line.startswith('>')]
            for line in contig_lines:
                contig_id = line.strip().replace('>', '')
                fasta_file_dict[index] = contig_id
                index += 1
            genome_id = fasta_file.split('.fasta')[0]
            fasta_dict[genome_id] = fasta_file_dict

    return fasta_dict


def get_contig_rows(contig, df, identical=False):
    """
    Given a contig and a dataframe containing contents of a BLAST file, returns a dataframe of all rows
    with that contig as either the query or sub.
    If identical, find rows where query and sub id are the same.
    """
    if identical:
        contig_rows = df.filter((df['query_id'] == contig) & (df['sub_id'] == contig))
    else:
        contig_rows = df.filter((df['query_id'] == contig) | (df['sub_id'] == contig))

    return contig_rows


def load_BLAST_file_df(file_path):
    """
    Loads a BLAST file assuming all data columns for output format 6 are present in default order.
    Ensures columns are loaded with correct variable types.
    """
    # Initialize empty dataframe to hold desired contig data
    df = pl.read_csv(file_path, has_header=False, sep='\t', new_columns=['query_id', 'sub_id', 'PI', 'bitscore'])
    df = df.with_columns(pl.col("query_id").cast(pl.Utf8, strict=False))
    df = df.with_columns(pl.col("sub_id").cast(pl.Utf8, strict=False))
    df = df.with_columns(pl.col("PI").cast(pl.Float64, strict=False))
    df = df.with_columns(pl.col("bitscore").cast(pl.Float64, strict=False))

    return df


def append_FASTA_neighborhood_contigs(filename, gene, contig_dict, identical=False):
    """
    Given a dictionary of contigs in a given neighborhood, where keys are int and values correspond to contigs,
    iterates through the contigs to obtain all relevant BLAST rows to a given dataframe.
    """
    # Dataframe to append to
    df = pl.DataFrame(schema={'query_id': pl.Utf8, 'sub_id': pl.Utf8, 'PI': pl.Float64, 'bitscore': pl.Float64})

    # Dataframe containing genome BLAST results
    data_df = load_BLAST_file_df(filename)
    # Concatenate all neighborhood row data using contig names
    for contig in contig_dict.keys():
        contig_rows = get_contig_rows(contig_dict[contig], data_df, identical=identical)
        df = pl.concat([df, contig_rows], how='vertical')

    del data_df

    return df


def make_assembly_dict(assembly_path):
    """
    Makes dictionary consisting of assembly file names with keys as the filename without the file extension, values
    as the full filename.
    """
    faa_dict = {}
    for filename in os.listdir(assembly_path):
        key = filename.split('.')[0]
        faa_dict[key] = filename

    return faa_dict


def get_blast_df(gene, blast_path, assembly_path, fasta_path, fasta_dict):
    """
    Loads whole genome BLAST into a dataframe with only relevant contig data rows.
    """
    # Identify all genomes the gene is present in using FASTA outputs
    faa_files = make_assembly_dict(assembly_path)
    present_genomes = [key for key in fasta_dict.keys()]

    blast_dict = dict()
    for genome_1, genome_2 in itertools.combinations(present_genomes, 2):

        # Initialize empty dataframe to hold desired contig data
        df = pl.DataFrame(schema={'query_id': pl.Utf8, 'sub_id': pl.Utf8, 'PI': pl.Float64, 'bitscore': pl.Float64})

        genome_1_filename = faa_files.get(genome_1)
        genome_2_filename = faa_files.get(genome_2)

        # Load comparative BLAST file
        if os.path.exists(blast_path + '/' + genome_1_filename + '.dmnd_' + genome_2_filename + '.txt'):
            file_pair_file_name = blast_path + '/' + genome_1_filename + '.dmnd_' + genome_2_filename + '.txt'
            contig_dict = fasta_dict[genome_1]
        elif os.path.exists(blast_path + '/' + genome_2_filename + '.dmnd_' + genome_1_filename + '.txt'):
            file_pair_file_name = blast_path + '/' + genome_2_filename + '.dmnd_' + genome_1_filename + '.txt'
            contig_dict = fasta_dict[genome_2]
        else:
            print("BLAST file for {g1} vs {g2} not found.".format(g1=genome_1_filename, g2=genome_2_filename))
            sys.exit(1)

        genome_1_file_name = blast_path + '/' + genome_1_filename + '.dmnd_' + genome_1_filename + '.txt'
        genome_2_file_name = blast_path + '/' + genome_2_filename + '.dmnd_' + genome_2_filename + '.txt'

        # Append all relevant BLAST rows for neighborhood contigs to df
        files = [file_pair_file_name, genome_1_file_name, genome_2_file_name]
        contig_dicts = [contig_dict, fasta_dict[genome_1], fasta_dict[genome_2]]
        identical_bools = [False, True, True]

        for blast_file_name, contig_dict, bool in zip(files, contig_dicts, identical_bools):
            blast_rows_df = append_FASTA_neighborhood_contigs(blast_file_name, gene, contig_dict,
                                                              identical=bool)
            df = pl.concat([df, blast_rows_df], how='vertical')

        blast_dict[genome_1 + '_' + genome_2 + '.blast.txt'] = df

    return blast_dict


def get_blast_dict_whole_genomes(assembly_path, fasta_path, output_path):
    """
    Creates dictionary for all extracted genes where each key is a gene name (str) and each value is the dictionary
    containing genome combinations that were blasted together as keys and their accompanying blast file data in a
    dataframe as values.
    """
    blast_path = '../../../' + output_path + '/diamond'
    extract_genes_list = os.listdir(fasta_path)

    BLAST_df_dict = dict()
    for gene in extract_genes_list:
        # Get the fasta file contents for each genome for that gene so contigs can be easily parsed
        fasta_dict = make_fasta_contig_dict(fasta_path, gene)

        # Make the BLAST dataframe using the contig dict
        occurrence_dict = get_blast_df(gene, blast_path, assembly_path, fasta_path, fasta_dict)

        BLAST_df_dict[gene] = occurrence_dict

    return BLAST_df_dict


def filter_BLAST_results(blast_gene_dict, cutoff=70):
    """
    Retains only protein sequences that share at least 70% overall sequence identity (i.e., homologs) in BLAST dict.
    """
    filtered_BLAST_dict = dict()
    for blast_file, blast_df in blast_gene_dict.items():
        # Filter rows to drop all rows with percent identity lower than specified threshold
        df = blast_df[blast_df['PI'] >= cutoff]

        # Remove rows where query ID and sub ID are the same: was only needed for normalized bitscores calculation
        final_df = df[df['query_id'] != df['sub_id']]

        filtered_BLAST_dict[blast_file] = df

    return filtered_BLAST_dict


def remove_BLAST_duplicates(blast_df):
    """
    Given a BLAST file dataframe, removes duplicate entries for same gene by keeping only the highest bitscore.
    """
    blast_df = blast_df.sort('normalized_bitscore', reverse=True)
    blast_df = blast_df.unique(subset=['query_id'], keep='first')

    return blast_df


def get_neighborhood_from_fasta(fasta_file_path):
    """
    Given a FASTA file of a gene neighborhood created using extraction.py, creates a list of all the gene
    identifiers in the neighborhood without their protein sequences.
    """
    neighborhood = []
    with open(fasta_file_path, 'r') as infile:
        contig_lines = [line for line in infile.readlines() if line.startswith('>')]
        for line in contig_lines:
            gene_id = line.replace('>', '').replace('\n', '')
            neighborhood.append(gene_id)

    return neighborhood


def get_neighborhoods_dict(fasta_dir_path):
    """
    Given user-provided path to FASTA files for the genomes being analyzed, returns a dictionary that contains the
    neighborhood of that gene in each genome for reference for later construction of similarity matrices.
    Keys: genes
    Values: Genome dict
           - Keys: Genome ID (str)
           - Values: Neighborhood
    """
    neighborhoods = dict()
    for dir in get_full_filepaths(fasta_dir_path):
        gene = get_filename(dir)
        gene_neighborhoods = dict()
        for fasta_filepath in get_full_filepaths(dir):
            fasta_filename = get_filename(fasta_filepath)
            neighborhood = get_neighborhood_from_fasta(fasta_filepath)
            gene_neighborhoods[fasta_filename] = neighborhood
        neighborhoods[gene] = gene_neighborhoods

    return neighborhoods


def get_similarity_matrices(blast_dict, neighborhoods_dict, neighborhood_size, fasta_path):
    """
    Calculates similarity matrix for each gene using data from BLAST dataframes and extracted neighborhoods
    """
    similarity_matrices_dict = dict()
    gene_genomes_dict = dict()

    for gene, blast_df_dict in blast_dict.items():

        # Get genome IDs for every genome that was used to generate the BLAST results
        genome_ids = [get_filename(filepath) for filepath in get_full_filepaths(fasta_path + '/' + gene)]

        # store in  gene genomes dict
        gene_genomes_dict[gene] = genome_ids

        # Set diagonal to max similarity score possible (for each neighborhood's score against itself)
        max_similarity_score = neighborhood_size * 2 + 1
        similarity_matrix = [[max_similarity_score for _ in range(len(genome_ids))] for _ in range(len(genome_ids))]

        for permutation in itertools.permutations(range(len(genome_ids)), 2):

            # Get genome indices and retrieve respective BLAST file data
            genome_1_index = permutation[0]
            genome_2_index = permutation[1]

            # Get genome names
            genome_1_id = genome_ids[genome_1_index]
            genome_2_id = genome_ids[genome_2_index]

            # Retrieve BLAST dataframe comparing the two genomes
            try:
                blast_filename = genome_1_id + '_' + genome_2_id + '.blast.txt'
                blast_df = blast_df_dict[blast_filename]

                # Get gene neighborhoods for comparison of neighborhood sizes (i.e. contig ends)
                neighborhood_1 = neighborhoods_dict[gene][genome_1_id]
                neighborhood_2 = neighborhoods_dict[gene][genome_2_id]

            except KeyError:
                try:
                    blast_filename = genome_2_id + '_' + genome_1_id + '.blast.txt'
                    blast_df = blast_df_dict[blast_filename]

                    # Get gene neighborhoods for comparison of neighborhood sizes (i.e. contig ends)
                    neighborhood_1 = neighborhoods_dict[gene][genome_2_id]
                    neighborhood_2 = neighborhoods_dict[gene][genome_1_id]
                except KeyError:
                    print("Key error for gene {g}: {a} and {b} didn't exist...".format(g=gene,
                                                                                       a=genome_1_id + '_' + \
                                                                                         genome_2_id + '.blast.txt',
                                                                                       b=genome_2_id + '_' + \
                                                                                         genome_1_id + '.blast.txt'))

            # ORIGINAL FILTERING CONDITION
            blast_neighborhood = blast_df.filter(((pl.col('query_id').is_in(neighborhood_1)) &
                                                  (pl.col('sub_id').is_in(neighborhood_2))) |
                                                 ((pl.col('query_id').is_in(neighborhood_2)) &
                                                  (pl.col('sub_id').is_in(neighborhood_1))))

            blast_neighborhood_df = remove_BLAST_duplicates(blast_neighborhood)

            sum_1 = round(blast_neighborhood_df['normalized_bitscore'].sum(), 3)

            if len(neighborhood_1) > len(neighborhood_2):
                len_diff = len(neighborhood_1) - len(blast_neighborhood_df)
                sum_2 = sum_1 + len(neighborhood_1) - len_diff - len(blast_neighborhood_df)

            elif len(neighborhood_2) > len(neighborhood_1):
                len_diff = len(neighborhood_2) - len(blast_neighborhood_df)
                sum_2 = sum_1 + len(neighborhood_2) - len_diff - len(blast_neighborhood_df)

            else:
                sum_2 = sum_1

            # Update values in upper and lower diagonal
            final_sum = round(sum_2, 3)
            similarity_matrix[genome_1_index][genome_2_index] = final_sum
            similarity_matrix[genome_2_index][genome_1_index] = final_sum

        similarity_matrices_dict[gene] = similarity_matrix

    return similarity_matrices_dict, gene_genomes_dict


def get_average_similarity_score(np_similarity_matrix, genome_names):
    """
    Computes average similarity score across all genomes for a given  gene's similarity matrix.
    """
    df = pd.DataFrame(data=np_similarity_matrix, index=genome_names, columns=genome_names)
    values = df.values

    return int(round(values.mean()))


def get_distance_matrix(np_similarity_matrix, genome_names):
    """"
    Converts numpy symmetric matrix array into a distance matrix dataframe!
    """
    df = pd.DataFrame(data=np_similarity_matrix, index=genome_names, columns=genome_names)
    values = df.values
    df = df.transform(lambda x: 1 - x / values.max())
    df = df.round(decimals=3)

    return df


def get_distance_std_dev(distance_matrix, genome_names=None):
    """
    Returns the standard deviation computed over all values across all rows and columns of the distance matrix.
    """
    df = pd.DataFrame(data=distance_matrix, index=genome_names, columns=genome_names)
    values = df.to_numpy()
    return np.std(values)


def get_minimum_distance_score(distance_matrix, genome_names=None):
    """
    Returns the smallest numerical value found in the distance matrix.
    """
    df = pd.DataFrame(data=distance_matrix, index=genome_names, columns=genome_names)
    values = df.to_numpy()
    return values.min()


def get_maximum_distance_score(distance_matrix, genome_names=None):
    """
    Returns the largest numerical value found in the distance matrix.
    """
    df = pd.DataFrame(data=distance_matrix, index=genome_names, columns=genome_names)
    values = df.to_numpy()
    return values.max()


def get_sparse_matrix(np_matrix):
    """
    Converts a numpy distance matrix to a scipy sparse matrix (important for use with MCL clustering).
    """
    sparse_matrix = sparse.csr_matrix(np_matrix)

    return sparse_matrix


def UPGMA_clustering(condensed_distance_matrix):
    """
    Applies UPGMA clustering to neighborhood similarity or symmetric distance matrix.
    """
    try:
        upgma_linkage = hierarchy.average(condensed_distance_matrix)
        return upgma_linkage
    except ValueError:
        print("Empty distance matrix was passed for the gene!")


def MCL_clustering(matrix, inflation):
    """
    Performs MCL clustering on a neighborhood sparse similarity or symmetric distance matrix.
    """
    sparse_matrix = get_sparse_matrix(matrix)
    result = mcl.run_mcl(sparse_matrix, inflation=inflation)
    clusters = mcl.get_clusters(result)

    return clusters


def DBSCAN_clustering(np_distance_matrix, epsilon, minpts):
    """
    Applies DBSCAN clustering to a neighborhood similarity or symmetric distance matrix.
    """
    distance_matrix = StandardScaler().fit_transform(np_distance_matrix)
    dbscan = DBSCAN(eps=epsilon, min_samples=minpts).fit(distance_matrix)

    return dbscan, dbscan.labels_


def check_clustering_savepaths(output_path):
    """
    Ensures that within the specified output directory, there is a directory for clustering visualizations with
    a distinct subdirectory containing the clustering results graph for each  gene with that algorithm.
    """
    for clustering_type in ['UPGMA', 'DBSCAN', 'MCL']:
        save_path = output_path + '/clustering/' + clustering_type
        check_output_path(save_path)


def get_num_clusters(clustering_labels):
    """
    Returns the number of unique clusters identified using a clustering algorithm given we provide a list of labels
    ('leaves_color_list' if upgma, clustering.labels_ if dbscan, clusters if mcl).
    """
    # If DBSCAN classifies all points as noise, we consider this equivalent to finding no clusters at all
    if -1 in set(clustering_labels) and len(set(clustering_labels)) == 1:
        return 0
    else:
        return len(set(clustering_labels))

def cluster_neighborhoods(assembly_path, fasta_path, blast_path, output_path,
                          neighborhood_size=10, inflation=2, epsilon=0.5, minpts=5):
    """
    Driver script for clustering neighborhoods obtained from extraction module.
    """

    # Make BLAST dataframe dictionary
    print("Fetching relevant BLAST data from DIAMOND outputs for each respective neighborhood...")
    BLAST_df_dict = get_blast_dict_whole_genomes(assembly_path, fasta_path, output_path)

    # Calculate normalized bitscores and save as column
    print("Calculating normalized bitscores for downstream scoring...")
    final_BLAST_dict = get_normalized_bitscores(BLAST_df_dict, blast_path)
    # final_BLAST_dict[gene] = filter_BLAST_results(blast_files_dict)

    # Update UID genes in JSON neighborhood representation
    # update_unidentified_genes_data(BLAST_df_dict, output_path, surrogates=False)

    # Update links in each  gene JSON file to reflect percent identities of blast hits
    print("Updating JSON neighborhood representations' links percent identities according to BLAST results...")
    update_JSON_links_PI(BLAST_df_dict, output_path, surrogates=False)
    update_JSON_links_PI(BLAST_df_dict, output_path, surrogates=True)

    # Get neighborhoods dict for calculating similarity matrices (needed to compare contig ends)
    for gene_subdir in os.listdir(fasta_path):
        # Remove non-FASTA files
        remove_files(fasta_path + '/' + gene_subdir, '.fasta')

    neighborhoods = get_neighborhoods_dict(fasta_path)

    # Calculate similarity matrix for each gene
    print("Calculating similarity matrices for each  gene...")
    check_output_path(output_path)
    similarity_matrices_dict, genome_names_dict = get_similarity_matrices(final_BLAST_dict,
                                                                          neighborhoods,
                                                                          neighborhood_size,
                                                                          fasta_path)
    # Make output dir for clustering results
    check_output_path(output_path + '/clustering')

    # Generate distance matrices
    print("Calculating distance matrices...")

    distance_matrices_df_dict = dict()
    average_similarity_scores_dict = dict()
    for gene, similarity_matrix_data in similarity_matrices_dict.items():

        similarity_matrix = np.array(similarity_matrix_data)
        genome_names = genome_names_dict[gene]

        try:
            distance_matrix = get_distance_matrix(similarity_matrix, genome_names)
            distance_matrices_df_dict[gene] = distance_matrix
            average_similarity_scores_dict[gene] = get_average_similarity_score(similarity_matrix, genome_names)

        # if directory empty of BLAST results
        except ValueError:
            print("Was unable to locate BLAST results for gene {g}. Skipped.".format(g=gene))
            pass

    print("Clustering neighborhoods...")
    check_clustering_savepaths(output_path)

    max_distance_scores_dict = dict()
    upgma_num_clusters_dict = dict()
    dbscan_num_clusters_dict = dict()
    mcl_num_clusters_dict = dict()

    os.mkdir(output_path + '/clustering/distance_matrices')
    os.mkdir(output_path + '/clustering/similarity_matrices')

    for gene, distance_matrix_df in distance_matrices_df_dict.items():

        # Get distance matrix
        np.fill_diagonal(distance_matrix_df.values, 0)
        distance_matrix = np.array(distance_matrix_df.values)

        genome_names = genome_names_dict[gene]

        # Retain maximum distance score for each  gene for distances histogram
        max_distance_scores_dict[gene] = get_maximum_distance_score(distance_matrix, genome_names)

        # Reorganize JSON clusters according to UPGMA order (left to right)
        genome_to_num_mapping = dict()
        for i in range(len(genome_names)):
            genome_to_num_mapping[str(i)] = genome_names[i]

        print("Generating UPGMA clusters for {g}...".format(g=gene))

        try:
            # Get UPGMA linkage matrix
            upgma_linkage = UPGMA_clustering(distance_matrix)

            # Use linkage matrix to build dendrogram used for updating JSON data: reorder clusters, make UPGMA view
            upgma_dendrogram = hierarchy.dendrogram(upgma_linkage)

            # Save number of clusters to upgma_num_clusters_dict
            num_clusters = get_num_clusters(upgma_dendrogram['leaves_color_list'])
            upgma_num_clusters_dict[gene] = num_clusters

            order_JSON_clusters_UPGMA(output_path, gene, upgma_dendrogram, genome_to_num_mapping, surrogates=False)
            order_JSON_clusters_UPGMA(output_path, gene, upgma_dendrogram, genome_to_num_mapping, surrogates=True)
            make_representative_UPGMA_cluster_JSON(output_path, gene, upgma_dendrogram, genome_to_num_mapping)

            # Interactive dendrogram visualization
            plotly_dendrogram(distance_matrix, genome_names, gene, output_path)

            # Save distance matrix as textfile
            check_output_path(output_path + '/clustering/distance_matrices')
            distance_matrix_df.to_csv(output_path + '/clustering/distance_matrices/' + gene + \
                                      '_distance_matrix.csv', sep='\t', index=False)

        except IndexError:
            print("Unable to perform UPGMA clustering for gene {g}. " \
                  "UPGMA results for {g} will be omitted.".format(g=gene))

        except TypeError:
            print("Unable to perform UPMA clustering for gene {g}. " \
                  "UPGMA results for {g} will be omitted.".format(g=gene))

        print("Generating DBSCAN clusters for {g}...".format(g=gene))

        try:
            dbscan_clusters, labels = DBSCAN_clustering(distance_matrix, epsilon, minpts)

            # Retain number of clusters
            num_dbscan_clusters = get_num_clusters(labels)
            dbscan_num_clusters_dict[gene] = num_dbscan_clusters

            # Plot DBSCAN clusters using distance matrix and PCoA: colour according to cluster assignment
            plotly_pcoa(distance_matrix_df, genome_names, labels, gene, output_path)

        except IndexError:
            print("Unable to perform DBSCAN clustering for gene {g}. " \
                  "DBSCAN results for {g} will be omitted.".format(g=gene))

        except KeyError:
            print("Unable to perform DBSCAN clustering for gene {g}. " \
                  "DBSCAN results for {g} will be omitted.".format(g=gene))

    for gene, similarity_matrix_data in similarity_matrices_dict.items():

        try:
            similarity_matrix = np.array(similarity_matrix_data)
            genome_names = genome_names_dict[gene]

            df = pd.DataFrame(data=similarity_matrix, index=genome_names, columns=genome_names)
            values = df.values

            print("Generating Markov clusters for {g}...".format(g=gene))

            clusters = MCL_clustering(similarity_matrix, inflation)

            # Retain number of clusters
            num_mcl_clusters = get_num_clusters(clusters)
            mcl_num_clusters_dict[gene] = num_mcl_clusters

            # Plot MCL clusters
            sparse_sim_matrix = get_sparse_matrix(similarity_matrix)
            plotly_mcl_network(sparse_sim_matrix, clusters, genome_names, gene, output_path)

            # Save distance matrix as textfile
            check_output_path(output_path + '/clustering/similarity_matrices')
            similarity_matrix_df = pd.DataFrame(data=similarity_matrix, index=genome_names, columns=genome_names)
            similarity_matrix_df.to_csv(output_path + '/clustering/similarity_matrices/' + gene + \
                                        '_similarity_matrix.csv', sep='\t', index=False)

        except IndexError:
            print("Unable to perform MCL clustering for gene {g}. " \
                  "MCL results for {g} will be omitted.".format(g=gene))

        except ValueError:
            print("Unable to perform MCL clustering for gene {g}. " \
                  "MCL results for {g} will be omitted.".format(g=gene))

    # Generate summary histograms for analyzed genomes' similarities
    print("Generating average similarity and max distance histograms for all neighborhoods...")
    plot_similarity_histogram(average_similarity_scores_dict, output_path)
    plot_distance_histogram(max_distance_scores_dict, output_path)
    make_clustering_summary_csv(output_path, distance_matrices_df_dict, upgma_num_clusters_dict,
                                mcl_num_clusters_dict, dbscan_num_clusters_dict)


def main(args=None):
    args = parse_args(args)
    cluster_neighborhoods(args.ASSEMBLY_PATH, args.FASTA_PATH, args.BLAST_PATH, args.OUTPUT_PATH,
                          args.n, args.i, args.e, args.m)


if __name__ == '__main__':
    sys.exit(main())