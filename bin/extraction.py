#!/usr/bin/env python

"""
Given extract files and their corresponding GBK files for complete bacterial genomes,
this script is used for identification of all unique gene neighborhoods present for easy cross-genome comparison.
"""

import os
import glob
import pandas as pd
from Bio import SeqIO
import json
import sys
import argparse
import itertools

from utils import get_filename, check_output_path, strip_brackets
from json_utils import (
    make_neighborhood_JSON_data,
    write_neighborhood_JSON,
    make_gene_HTML,
)


def parse_args(args=None):
    Description = "Extract AMR gene neighborhoods according to fixed window size of N genes upstream and downstream."
    Epilog = (
        "Example usage: python extraction.py <extract_PATH> <GBK_PATH> <OUTPUT_PATH> -n <NEIGHBORHOOD_SIZE> -p "
        "<PERCENT> "
    )

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-i",
        dest="INPUT_FILE_PATH",
        metavar="in_file",
        type=str,
        help="Path to an input textfile specifying "
        "either a) column names for data needed for"
        "extracting from pre-existing annotation "
        "files or b) gene names to extract.",
    )
    parser.add_argument(
        "-x",
        dest="EXTRACT_PATH",
        metavar="extract",
        type=str,
        help="Path to directory containing extract "
        "files. Must have names corresponding to "
        "GBK files.",
        nargs="*",
    )
    parser.add_argument(
        "-g",
        dest="GBK_PATH",
        metavar="gbk",
        type=str,
        help="Path to directory containing GBK files. \
                                                                   Must have names corresponding to extract files.",
        nargs="*",
    )
    parser.add_argument(
        "-o",
        dest="OUTPUT_PATH",
        metavar="output_path",
        type=str,
        help="Path to output directory where "
        "extracted neighborhood FASTA files will"
        " be saved.",
    )
    parser.add_argument(
        "-w",
        dest="HTML_TEMPLATE",
        metavar="html_template",
        type=str,
        help="Path to HTML template.",
    )
    parser.add_argument(
        "-n",
        metavar="n",
        type=int,
        default=10,
        help="Neighborhood window size, i.e. number of genes "
        "to consider upstream and downstream of focal "
        "gene.",
    )
    parser.add_argument(
        "-p",
        metavar="p",
        type=float,
        default=0.75,
        help="Cutoff percentage of genomes that a given "
        "AMR gene should be present in for its "
        "neighborhood to be considered.",
    )
    parser.add_argument(
        "-c",
        metavar="label_cols",
        type=str,
        default=None,
        help="If using annotation files predicting "
        "features, list of space separated "
        "column names to be added to the gene"
        " names.",
    )
    return parser.parse_args(args)


def load_input_file(input_file_path):
    """
    Handles user input file.

    I) If a list of column specifiers for custom annotation files in format:

    annotation_contig_column_name = Contig
    annotation_gene_name_column_name = Gene_Name
    annotation_gene_start_column_name = Gene_Start
    annotation_gene_end_column_name = Gene_End

    Returns a dictionary where keys are the column names, values are the equivalent columns in the custom
    annotation files.

    II) If a list of gene names in format:

    gene_A
    gene_B
    ...
    gene_X

    Returns a list of non-redundant gene names to extract.
    """
    with open(input_file_path, "r") as infile:
        input_data = infile.readlines()
        if len(input_data) == 0:
            print(
                "Error: No data found in input file. Please provide an input file specifying either a) column names "
                "of required data columns in custom annotation files, or b) list of genes to extract, in specified "
                "format given in docs."
            )
            sys.exit(1)

    # Determine which input file format was given
    if "=" in input_data[0]:
        # Custom annotation files provided
        col_names_dict = dict()
        for line in input_data:
            col_name_tokens = line.split("=")
            col_names_dict[col_name_tokens[0].strip(" ")] = (
                col_name_tokens[1].strip(" ").strip("\n")
            )
        return col_names_dict

    else:
        # List of gene names provided
        gene_names = []
        for line in input_data:
            gene_names.append(line.strip("\n"))
        return gene_names


def get_gbk_data(gbk_filepaths):
    """
    Return dictionary of dataframes containing each GBK file's Genbank records, features and annotations.
    """
    # Contains a dataframe for each GBK file, where keys are filenames
    gbk_dataframes = {}

    for gbk_filepath in gbk_filepaths:
        # Get GBK filename for GBK dictionary keys
        gbk_filename = get_filename(gbk_filepath)

        # Get GBK dataframe and contig names
        gbk_df = make_GBK_dataframe(gbk_filepath)

        # Preprocess all bracketed columns to remove brackets
        for col in gbk_df:
            if type(gbk_df[col][0]) is list or isinstance(gbk_df[col][0], str):
                gbk_df[col] = gbk_df[col].apply(lambda x: strip_brackets(x))
                gbk_df[col] = gbk_df[col].apply(lambda x: x.strip("'"))

        # Retain GBK dataframe
        gbk_dataframes[gbk_filename] = gbk_df

    return gbk_dataframes


def load_GBK_file(GBK_filepath):
    """
    Loads a GBK file and gets sequence records, features, and annotations for the file
    """
    filename, extension = os.path.splitext(GBK_filepath)
    assert (
        extension == ".gbk" or extension == ".gb"
    ), "Error: filepath provided does not lead to a Genbank file."

    # Extract and store all features and annotations per record
    records = [record for record in SeqIO.parse(GBK_filepath, "genbank")]
    features = [record.features for record in records]
    annotations = [record.annotations for record in records]

    return records, features, annotations


def make_GBK_dataframe(GBK_file_path):
    """
    Input: The output file from parse_genbank_proteins
    Returns a dataframe containing GBK data
    """
    gbk_df = pd.DataFrame()

    gene_start = []
    gene_end = []
    gene_strand = []
    gene_name = []
    loc_tag = []
    function = []
    protein_seq = []
    contig_name = []

    for index, record in enumerate(SeqIO.parse(GBK_file_path, "genbank")):
        for feature in record.features:
            if feature.type == "CDS" and "pseudo" not in feature.qualifiers:
                gene_start.append(feature.location.start)
                gene_end.append(feature.location.end)
                gene_strand.append(feature.location.strand)
                loc_tag.append(feature.qualifiers["locus_tag"])
                function.append(feature.qualifiers["product"])
                protein_seq.append(str(feature.qualifiers["translation"]))
                contig_name.append(record.id)

                # Preserve gene name if present, otherwise mark as UID (unidentified)
                if "gene" in feature.qualifiers:
                    gene_name.append(feature.qualifiers["gene"])
                else:
                    gene_name.append(
                        "UID-"
                        + str(feature.qualifiers["locus_tag"]).strip("[").strip("]")
                    )

    gbk_df["Gene_Start"] = gene_start
    gbk_df["Gene_End"] = gene_end
    gbk_df["Gene_Strand"] = gene_strand
    gbk_df["Locus_Tag"] = loc_tag
    gbk_df["Gene_Name"] = gene_name
    gbk_df["Product"] = function
    gbk_df["Protein_Sequence"] = protein_seq
    gbk_df["Contig_Name"] = contig_name

    return gbk_df


def make_extraction_dataframe(input_file_path, col_names_data):
    """
    Input: Path to tab-delimited file listing custom annotation data.
    Returns a dataframe containing
    """
    extraction_df = pd.read_csv(input_file_path, sep="\t", header=0)

    # Drop unnecessary columns
    drop_cols = [
        col_name
        for col_name in extraction_df.columns.values.tolist()
        if col_name not in col_names_data.keys()
    ]

    # Rename columns
    extraction_df.rename(columns=col_names_data, inplace=True)

    return extraction_df


def clean_gene_identifier(name):
    """
    Simplifies gene names and removes extraneous characters. Mostly suitable for RGI ARO name cleanup currently.
    (e.g. 'mecC-type BlaZ' -> 'mecC_BlaZ',
          'vanS gene in vanN cluster' -> 'vanS',
          'Bifidobacterium adolescentis rpoB mutants conferring resistance to rifampicin' -> 'rpoB')
    """
    tokenized_name = name.split(" ")
    try:
        if len(tokenized_name) > 1:
            short_name = (
                tokenized_name[0][0] + tokenized_name[1][0] + "_" + tokenized_name[2]
            )
        else:
            short_name = tokenized_name[0]
    except IndexError:
        short_name = tokenized_name[0]

    short_name = (
        short_name.replace("'", "")
        .replace("-", "_")
        .replace("(", "")
        .replace(")", "")
        .replace("/", "-")
        .replace(" ", "_")
        .replace(".", "-")
        .replace('"', "")
        .strip("-")
        .strip('"')
    )

    return short_name


def clean_gene_names_df(extract_df, genome_name):
    """
    Simplifies gene names from annotation df
    """
    extract_df.reset_index(drop=True, inplace=True)
    original_names = extract_df["Gene_Name"].tolist()
    extract_df["Gene_Name"] = extract_df["Gene_Name"].apply(
        lambda name: clean_gene_identifier(name)
    )
    final_names = extract_df["Gene_Name"].tolist()
    names_dict = dict(zip(original_names, final_names))
    return extract_df, names_dict


def partition_contig_name(contig_str):
    """
    Splits contig name string and retains first part only (i.e. contig number)
    """
    return contig_str.split("_")[0]


def process_GBK_df_for_BLAST(gbk_df):
    """
    Given a dataframe for a GBK file, removes bracket formatting for easier processing and readability later
    """
    gbk_df["Locus_Tag"] = gbk_df["Locus_Tag"].apply(
        lambda i: str(i).replace("[", "").replace("]", "").replace("'", "")
    )
    gbk_df["Gene_Name"] = gbk_df["Gene_Name"].apply(
        lambda i: str(i)
        .replace("[", "")
        .replace("]", "")
        .replace("'", "")
        .replace('"', "")
        .replace("/", "")
        .replace("(", "")
        .replace(")", "")
    )
    gbk_df["ProteinSequence"] = gbk_df["ProteinSequence"].str.strip("[]")

    return gbk_df


def make_extract_df_contig_col(extract_df):
    """
    Applies transformation to contig name column to retain only first part of default contig col value from GBK.
    """
    extract_df = extract_df["Contig"].apply(
        lambda contig_name: contig_name.split("_")[0]
    )

    return extract_df


def group_by_contig(gbk_df):
    """
    Creates a dict where keys are contig identifiers, and values correspond to all GBK data present for that contig
    """
    group = gbk_df.groupby(gbk_df["Contig_Name"])
    datasets = {}

    for groups, data in group:
        datasets[groups] = data

    return datasets


def get_unique_AMR_genes(extract_hit_dicts):
    """
    Helper function to return a list of all unique AMR genes present across all inputted genomes.
    """
    unique_AMR_genes = set()
    for extract_hit_dict in extract_hit_dicts:
        # Get a list of all AMR genes present in the genome
        AMR_genes = list(extract_hit_dict.keys())
        # Add to set
        unique_AMR_genes.add(tuple(AMR_genes))

    return unique_AMR_genes


def find_union_genes(extract_dataframes):
    """
    Returns a list of unique genes occurring in the full set of genomes from the extract dataframes.
    """
    unique_genes = []
    for genome_id, extract_df in extract_dataframes.items():
        genome_genes = extract_df["Gene_Name"].to_list()
        for gene in genome_genes:
            unique_genes.append(gene)

    return list(set(unique_genes))


def make_occurrence_dict(extract_dataframes, gene):
    """
    Given the extract dataframes for the genomes being analyzed, creates a dictionary of AMR genes
    """
    gene_occurrences_dict = {}
    for genome, df in extract_dataframes.items():
        gene_row = df.loc[df["Gene_Name"] == gene]
        if len(gene_row) > 0:
            gene_occurrences_dict[genome] = gene_row

    return gene_occurrences_dict


def make_gene_neighborhood_df(
    GBK_df_dict, genome_id, gene_start, gene_name, neighborhood_size, modified_gene_name
):
    """
    Finds gene neighborhood of size 2N (user-defined by neighborhood_size, i.e., N genes downstream and upstream)
    for a given AMR gene for cross genome comparison.
    """
    # Get the GBK data for the given genome
    try:
        GBK_df = GBK_df_dict[genome_id]
        GBK_df.reset_index(drop=True, inplace=True)
    except KeyError:
        return TypeError

    # Keep track of contig ends: track relative start and stop coordinates of neighborhood
    neighborhood_indices = []

    # Get the focal (central) AMR gene
    # Subtract one from gene start index to account for automatic padding
    gene_df_row = GBK_df.loc[((GBK_df["Gene_Start"] == gene_start - 1))]
    if gene_df_row.empty:
        gene_df_row = GBK_df.loc[((GBK_df["Gene_Start"] == gene_start))]

    try:
        # Get only genes on the same contig as the focal gene for consideration as neighbors
        contig_id = gene_df_row.Contig_Name.tolist()[0]
        contig_df = GBK_df.loc[(GBK_df["Contig_Name"] == contig_id)].sort_values(
            by="Gene_Start"
        )
        contig_df.reset_index(drop=True, inplace=True)
    except IndexError:
        return TypeError

    gene_index_coord = contig_df.index[
        (contig_df["Gene_Start"] == gene_start - 1)
    ].tolist()
    if not gene_index_coord:
        gene_index_coord = contig_df.index[
            (contig_df["Gene_Start"] == gene_start)
        ].tolist()

    gene_index = gene_index_coord[0]

    gene_df_row = contig_df.iloc[[gene_index]]
    if modified_gene_name is not None:
        gene_df_row.loc[:, "Gene_Name"] = modified_gene_name

    # Get downstream neighbors
    downstream = [gene_index - index for index in range(1, neighborhood_size + 1)]

    # If contig end present downstream, some indices will be negative: remove these to prevent index errors
    downstream_indices = [index for index in downstream if index >= 0]
    downstream_neighbors = pd.DataFrame(
        columns=[
            "Gene_Start",
            "Gene_End",
            "Gene_Strand",
            "Locus_Tag",
            "Gene_Name",
            "Product",
            "Protein_Sequence",
            "Contig_Name",
        ]
    )
    for i in range(len(downstream_indices) - 1, -1, -1):
        try:
            neighbor = contig_df.iloc[[downstream_indices[i]]]
            downstream_neighbors = pd.concat([downstream_neighbors, neighbor])
        except IndexError:
            print("Contig end found at position -{} downstream.".format(i + 1))
            neighborhood_indices.append(i)

    # If there was no contig end, append default N size to neighborhood_indices
    if len(neighborhood_indices) == 0:
        neighborhood_indices.append(-neighborhood_size)

    # Get upstream neighbors
    upstream_indices = [gene_index + index for index in range(1, neighborhood_size + 1)]
    upstream_neighbors = pd.DataFrame(
        columns=[
            "Gene_Start",
            "Gene_End",
            "Gene_Strand",
            "Locus_Tag",
            "Gene_Name",
            "Product",
            "Protein_Sequence",
            "Contig_Name",
        ]
    )
    for i in range(len(upstream_indices)):
        contig_end_found = False
        try:
            neighbor = contig_df.iloc[[upstream_indices[i]]]
            upstream_neighbors = pd.concat([upstream_neighbors, neighbor])
        except IndexError:
            if not contig_end_found:
                print("Contig end found at position {} upstream.".format(i + 1))
                contig_end_found = True
            if len(neighborhood_indices) < 2:
                neighborhood_indices.append(i)

    if len(neighborhood_indices) < 2:
        neighborhood_indices.append(neighborhood_size)

    neighborhood_df = pd.concat([downstream_neighbors, gene_df_row, upstream_neighbors])
    return neighborhood_df, neighborhood_indices


def get_all_gene_neighborhoods(
    gene_instance_dict, GBK_df_dict, unique_genes, neighborhood_size, cols
):
    """
    Given a dictionary of genes determined using extract outputs and a list of unique genes, creates one dictionary
    containing all gene neighborhoods for a fixed size of 2N (user-defined by neighborhood_size, i.e., N genes
    downstream and upstream) for the set of genomes being analyzed.
    """
    # Will be used to store dataframes of neighborhoods
    gene_neighborhoods = {}

    # Keeps track of contig ends, if applicable
    contig_ends = {}

    for gene, gene_dict in gene_instance_dict.items():
        # Keep track of size 2N neighborhood
        neighbors = {}

        # Keep track of start and stop indices for each neighborhood
        contig_end_flags = {}
        errors = 0

        for genome, data in gene_dict.items():
            # Get start index of the focal AMR gene from the extract dataframe
            start_vals = gene_dict[genome]["Gene_Start"]
            start_vals_list = list(start_vals)
            start_index = start_vals_list[0]

            # If using annotation, set gene name according to annotated feature name
            if cols is not None:
                additional_data = []
                gene_dict[genome].reset_index(drop=True, inplace=True)
                for col in cols:
                    additional_data.append(gene_dict[genome].loc[0, col])
                delim = "_"
                temp_data_str = list(map(str, additional_data))
                str_data = delim.join(temp_data_str)
                modified_gene_name = str(gene + "_" + str_data)
            else:
                modified_gene_name = None

            # Make gene neighborhood dataframe for each genome for the focal gene, AMR_gene
            try:
                neighborhood_df, indices = make_gene_neighborhood_df(
                    GBK_df_dict,
                    genome,
                    start_index,
                    gene,
                    neighborhood_size,
                    modified_gene_name,
                )
                neighborhood_df.reset_index(drop=True, inplace=True)

                if len(neighborhood_df) > 1:
                    neighbors[genome] = neighborhood_df
                    contig_end_flags[genome] = indices
            except TypeError:
                errors += 1
                print("Gene {} was not present in this genome!".format(gene))

        gene_neighborhoods[gene] = neighbors
        contig_ends[gene] = contig_end_flags

    return gene_neighborhoods, contig_ends


def get_neighborhood_gene_data(neighborhood_df):
    """
    Given a neighborhood dataframe, obtains the locus tags, protein sequences, and gene names as separate dicts
    """
    locus_list = neighborhood_df["Locus_Tag"].tolist()
    protein_list = neighborhood_df["Protein_Sequence"].tolist()
    gene_name_list = neighborhood_df["Gene_Name"].tolist()

    locus_to_protein_dict = {}
    for gene in range(len(locus_list)):
        locus_to_protein_dict[locus_list[gene].strip()] = protein_list[gene]

    return locus_list, protein_list, gene_name_list


def get_neighborhood_data(neighborhoods_dict, num_neighbors):
    """
    Extracts locus, protein sequence, and gene name data for a dictionary of neighborhood data created using
    get_all_gene_neighborhoods.
    """
    locus_dict = {}
    protein_dict = {}
    gene_name_dict = {}

    for AMR_gene, genome_neighborhoods in neighborhoods_dict.items():
        locus_tags = {}
        protein_seqs = {}
        gene_names = {}

        for genome_id, neighborhood_df in genome_neighborhoods.items():
            (
                locus_tags[genome_id],
                protein_seqs[genome_id],
                gene_names[genome_id],
            ) = get_neighborhood_gene_data(genome_neighborhoods[genome_id])

        locus_dict[AMR_gene] = locus_tags
        protein_dict[AMR_gene] = protein_seqs
        gene_name_dict[AMR_gene] = gene_names

    return locus_dict, protein_dict, gene_name_dict


def write_gene_neighborhood_to_FNA(
    AMR_gene_neighborhoods_dict, AMR_gene, locuses_dict, protein_seqs_dict, out_path
):
    """
    Given a dictionary containing AMR gene neighborhoods for a set of genomes being analyzed and a specified output
    directory path, creates a distinct .fna file for each AMR gene neighborhood to use for BLAST All-vs-All comparison.
    """
    for genome_id, genome_neighborhood_df in AMR_gene_neighborhoods_dict.items():
        # Create a new FASTA file to write neighborhood sequence and identifier data to
        with open(out_path + "/" + genome_id + ".fasta", "w+") as output_fasta:
            for locus, protein_seq in zip(
                locuses_dict[AMR_gene][genome_id],
                protein_seqs_dict[AMR_gene][genome_id],
            ):
                locus_output = locus.strip("'")
                protein_output = protein_seq.strip("'")
                output_fasta.write(">{}\n".format(locus_output))
                output_fasta.write("{}\n".format(protein_output))
        output_fasta.close()


def delete_low_occurring_genes(
    gene_occurrence_dict, num_genomes, cutoff_percentage=0.25
):
    """
    Removes keys of genes not present in cutoff percentage of genomes.
    Cutoff percentage is a float between 0 and 1 (e.g., 0.25 means only genes present in min 25% of genomes are kept).
    """
    # Define threshold amount of genomes an AMR gene must be present in for inclusion
    minimum_num_genomes = round(num_genomes * cutoff_percentage)

    # Remove all AMR genes that occur in less genomes than the threshold amount
    genes_to_remove = []
    for gene, neighborhoods_dict in gene_occurrence_dict.items():
        if len(neighborhoods_dict) < minimum_num_genomes:
            genes_to_remove.append(gene)

    for gene in genes_to_remove:
        del gene_occurrence_dict[gene]

    return gene_occurrence_dict


def get_AMR_gene_statistics(extract_dataframes):
    """
    Given RGI dataframes, counts the number of Perfect and Strict hits in each genome.
    Needed for extraction summary file.
    """
    amr_gene_statistics = {}
    for genome, extract_df in extract_dataframes.items():
        try:
            perfect_count = len(extract_df[extract_df["Cut_Off"] == "Perfect"])
            strict_count = len(extract_df[extract_df["Cut_Off"] == "Strict"])
            loose_count = len(extract_df[extract_df["Cut_Off"] == "Loose"])
            amr_gene_statistics[genome] = [perfect_count, strict_count, loose_count]
        except IndexError:
            pass

    return amr_gene_statistics


def write_summary_file(
    output_dir,
    gbk_dataframes,
    neighborhood_size,
    extraction_genes,
    gene_name_identifiers_dict=None,
):
    """
    Writes summary data of neighborhood extraction for the set of genomes being analyzed to
    output/Neighborhoods_Extraction_Summary.txt.
    """
    with open(
        output_dir + "/" + "Neighborhoods_Extraction_Summary.txt", "w"
    ) as outfile:
        outfile.write(
            "----------------------------------------------------------------------------------------\n"
        )
        outfile.write("NEIGHBORHOOD EXTRACTION SUMMARY" + "\n")
        outfile.write(
            "----------------------------------------------------------------------------------------"
            + "\n\n"
        )
        outfile.write(
            "----------------------------------------------------------------------------------------\n"
        )
        outfile.write("General Details" + "\n")
        outfile.write(
            "----------------------------------------------------------------------------------------\n"
        )
        outfile.write(
            "Number of genomes analyzed: {}".format(len(gbk_dataframes.keys())) + "\n"
        )
        outfile.write(
            "Neighborhood size extracted: {}".format(neighborhood_size) + "\n"
        )
        outfile.write(
            "Number of genes whose neighborhoods were extracted: {}\n\n".format(
                len(extraction_genes)
            )
        )
        if gene_name_identifiers_dict is not None:
            outfile.write(
                "----------------------------------------------------------------------------------------"
                + "\n"
            )
            outfile.write("Original Gene Names and Simplified Names" + "\n")
            outfile.write(
                "----------------------------------------------------------------------------------------"
                + "\n"
            )
            for gene_name, simplified_gene_name in sorted(
                gene_name_identifiers_dict.items()
            ):
                outfile.write(gene_name + "\t--->\t" + simplified_gene_name + "\n")


def extract_neighborhoods(
    input_file_path,
    extract_path,
    gbk_path,
    output_path,
    html_template,
    num_neighbors,
    cutoff_percent,
    label_cols=None,
):
    """
    Driver script for extracting all gene neighborhoods from either a) specified GBK and extract files or b) GBK files
    according to provided list of gene names to output FASTA format protein sequences for each neighborhood.
    """
    # Load input file: determine case
    input_file_data = load_input_file(input_file_path)

    # Make output directory if non-existent
    check_output_path(output_path)

    # Case 1: User provided custom annotation files to extract data according to
    if type(input_file_data) is dict:
        # 0) Get user path args from pipeline for path to extract files and GBK files
        try:
            extract_filepaths, gbk_filepaths = extract_path, gbk_path
        except IndexError:
            print(
                "Please provide paths to the extract files and GBK files respectively when running this script."
            )
            sys.exit(1)

        # 1) Get Genbank records, features and annotations for each GBK file
        print("Processing GBK files...")
        gbk_dataframes = get_gbk_data(gbk_filepaths)

        # 2) Get extract dataframes and do all required preprocessing for later manipulation and comparison
        # Contains a dataframe for each extract file, where keys are filenames
        print("Processing extract files...")
        extract_dataframes = {}
        for extract_filepath in extract_filepaths:
            extract_filename = get_filename(extract_filepath, extract=True)
            extract_df = make_extraction_dataframe(extract_filepath, input_file_data)
            extract_dataframes[extract_filename] = extract_df

        gene_name_identifiers_dict = dict()
        for extract_filename, extract_df in extract_dataframes.items():
            # Add locus tag and contig details to the extract dataframes
            extract_dataframes[extract_filename] = make_extract_df_contig_col(
                extract_df
            )
            extract_dataframes[extract_filename], names_dict = clean_gene_names_df(
                extract_df, extract_filename
            )
            gene_name_identifiers_dict.update(names_dict)

            # Preprocess all bracketed columns to remove brackets
            for col in extract_df:
                if type(extract_df[col][0]) is list or isinstance(
                    extract_df[col][0], str
                ):
                    try:
                        extract_df[col] = extract_df[col].apply(
                            lambda x: strip_brackets(x)
                        )
                    except AttributeError:
                        pass

        # 4) Get unique gene instances for genes to extract
        extraction_genes = find_union_genes(extract_dataframes)

        # 5) Get occurrence dict for all genomes
        genomes_occurrence_dict = {}
        for gene in extraction_genes:
            # Make entry for gene occurrences
            genomes_occurrence_dict[gene] = make_occurrence_dict(
                extract_dataframes, gene
            )

        genomes_occurrence_dict_filtered = delete_low_occurring_genes(
            genomes_occurrence_dict, len(gbk_filepaths), cutoff_percent
        )

        print(
            "Extracting gene neighborhood data for neighborhoods of size {}...".format(
                num_neighbors
            )
        )

        # Process additional cols from annotations to add to gene name
        if label_cols is not None:
            col_tokens = label_cols.split(",")
            for col_token in col_tokens:
                col_token = col_token.replace(" ", "")
        else:
            col_tokens = None

        # 6) Extract gene neighborhoods and store them in dataframes
        neighborhoods, neighborhood_indices = get_all_gene_neighborhoods(
            genomes_occurrence_dict_filtered,
            gbk_dataframes,
            extraction_genes,
            num_neighbors,
            col_tokens,
        )

        # 7) Get the locus, protein, and gene name details for each neighborhood respectively for FNA file creation
        locuses, protein_seqs, gene_names = get_neighborhood_data(
            neighborhoods, num_neighbors
        )

        # 8) Save neighborhoods to FNA files needed for All-vs-all BLAST results later
        print("Generating neighborhood FNA files...")
        for gene, genome_neighborhoods in neighborhoods.items():
            # Make output subdirectory for the gene
            gene_name = (
                gene.replace("(", "_")
                .replace(")", "_")
                .replace("'", "")
                .replace("/", "")
            )
            out_path = output_path + "/fasta/" + gene_name

            check_output_path(out_path)

            for genome_id in genome_neighborhoods.keys():
                # Make file with genome ID as filename
                write_gene_neighborhood_to_FNA(
                    genome_neighborhoods, gene, locuses, protein_seqs, out_path
                )

        # 9) Make neighborhoods summary textfile in output dir
        print("Making extraction summary file...")
        write_summary_file(
            output_path,
            gbk_dataframes,
            num_neighbors,
            extraction_genes,
            gene_name_identifiers_dict,
        )

        # 10) Save gene neighborhoods and indices in textfile: needed for JSON representations downstream
        print("Generating neighborhood JSON representations...")
        for AMR_gene, neighborhood_data in neighborhoods.items():
            # Get data needed to write JSON files
            neighborhood_JSON_dict = make_neighborhood_JSON_data(
                neighborhood_data, AMR_gene, num_neighbors
            )

            # Create JSON file
            write_neighborhood_JSON(neighborhood_JSON_dict, AMR_gene, output_path)

        make_gene_HTML(neighborhoods.keys(), html_template, output_path)

        with open(output_path + "/" + "neighborhood_indices.json", "w+") as outfile:
            outfile.write(json.dumps(neighborhood_indices, indent=4, sort_keys=True))

        print("Neighborhood extraction complete.")

    # Case 2: User provided list of genes to extract
    elif type(input_file_data) is list:
        try:
            gbk_filepaths = glob.glob(os.path.join(gbk_path, "*.gbk"))
        except FileNotFoundError:
            print(
                "Error: there are no GBK files found in the specified directory. "
                "Please double-check the provided path."
            )
            sys.exit(1)

        # 1) Get Genbank records, features and annotations for each GBK file
        print("Processing GBK files...")
        gbk_dataframes = get_gbk_data(gbk_filepaths)

        # 2) Determine all occurrences of given genes in the genome dataset
        genomes_occurrence_dict = {}
        for gene in input_file_data:
            # Make entry for gene occurrences
            genomes_occurrence_dict[gene] = make_occurrence_dict(gbk_dataframes, gene)

        # 3) Remove low occurrence genes according to cutoff threshold
        genomes_occurrence_dict_filtered = delete_low_occurring_genes(
            genomes_occurrence_dict, len(gbk_filepaths), cutoff_percent
        )

        print(
            "Extracting gene neighborhood data for neighborhoods of size {}...".format(
                num_neighbors
            )
        )

        # 6) Extract gene neighborhoods and store them in dataframes
        neighborhoods, neighborhood_indices = get_all_gene_neighborhoods(
            genomes_occurrence_dict_filtered,
            gbk_dataframes,
            input_file_data,
            num_neighbors,
            label_cols,
        )

        # 7) Get the locus, protein, and gene name details for each neighborhood respectively for FNA file creation
        locuses, protein_seqs, gene_names = get_neighborhood_data(
            neighborhoods, num_neighbors
        )

        # 8) Save neighborhoods to FNA files needed for All-vs-all BLAST results later
        print("Generating neighborhood FNA files...")
        for gene, genome_neighborhoods in neighborhoods.items():
            # Make output subdirectory for the gene
            gene_name = gene.replace("(", "_").replace(")", "_").replace("'", "")
            out_path = output_path + "/fasta/" + gene_name
            check_output_path(out_path)

            for genome_id in genome_neighborhoods.keys():
                # Make file with genome ID as filename
                write_gene_neighborhood_to_FNA(
                    genome_neighborhoods, gene, locuses, protein_seqs, out_path
                )

        # 9) Make neighborhoods summary textfile in output dir
        print("Making extraction summary file...")
        write_summary_file(output_path, gbk_dataframes, num_neighbors, input_file_data)

        # 10) Save gene neighborhoods and indices in textfile: needed for JSON representations downstream
        print("Generating neighborhood JSON representations...")
        for gene, neighborhood_data in neighborhoods.items():
            # Get data needed to write JSON files
            neighborhood_JSON_dict = make_neighborhood_JSON_data(
                neighborhood_data, gene, num_neighbors
            )

            # Create JSON file
            write_neighborhood_JSON(neighborhood_JSON_dict, gene, output_path)

        make_gene_HTML(neighborhoods.keys(), html_template, output_path)

        with open(output_path + "/" + "neighborhood_indices.txt", "w+") as outfile:
            outfile.write(str(neighborhood_indices))

        print("Neighborhood extraction complete.")

    else:
        print(
            "Invalid input file provided. Please double-check the documentation for the proper input file format."
        )
        sys.exit(1)


def main(args=None):
    args = parse_args(args)
    extract_neighborhoods(
        args.INPUT_FILE_PATH,
        args.EXTRACT_PATH,
        args.GBK_PATH,
        args.OUTPUT_PATH,
        args.HTML_TEMPLATE,
        args.n,
        args.p,
        args.c,
    )


if __name__ == "__main__":
    sys.exit(main())
