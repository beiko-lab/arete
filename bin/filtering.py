#!/usr/bin/env python

"""
Given RGI files and their corresponding GBK files for complete bacterial genomes, as well as gene neighborhood clusters
generated using neighborhood_clustering.py, the utilities in this class are used for two main types of filtering:

a) Filtering of redundant AMR gene neighborhoods based on the clustering results for cleaner visualization downstream.
b) If user opts to investigate Loose hits using predictive module, filtering of low percent identity Loose hits.
"""

import os
import json
import polars as pl
from utils import check_output_path
from utils import make_fasta_contig_dict
from json_utils import load_JSON_data, get_focal_index, make_assembly_dict


def filter_genes_percent_identity(
    output_path, fasta_path, extract_gene, blast_df_gene_dict, num_neighbors
):
    """
    Used within clustering module in order to go over previous filtered neighborhoods representation and remove
    identical neighborhoods based on percent identity (PI) of matched hits.
    """
    # Load filtered neighborhood
    json_data, gene_path = load_JSON_data(output_path, extract_gene)

    # Go over combinations of genomes
    genome_genes_dict = dict()
    for cluster in json_data["clusters"]:
        genome_genes_dict[cluster["name"]] = False

    representative_genomes = dict()
    candidate_rep_genomes = [cluster["name"] for cluster in json_data["clusters"]]
    genomes = [cluster["name"] for cluster in json_data["clusters"]]

    for genome in candidate_rep_genomes:
        represented_genomes = []
        for genome_2 in genomes:
            if (
                genome_2 not in representative_genomes.keys()
                and genome_genes_dict[genome_2] == False
            ):
                equivalent = True
                if genome != genome_2:
                    cluster_1 = [
                        cluster
                        for cluster in json_data["clusters"]
                        if cluster["name"] == genome
                    ][0]
                    cluster_2 = [
                        cluster
                        for cluster in json_data["clusters"]
                        if cluster["name"] == genome_2
                    ][0]
                    clust_1 = cluster_1["loci"][0]
                    clust_2 = cluster_2["loci"][0]

                    # Get genome IDs
                    genome_1_id = cluster_1["name"]
                    genome_2_id = cluster_2["name"]

                    if (
                        genome_1_id + "_" + genome_2_id + ".blast.txt"
                        in blast_df_gene_dict
                    ):
                        blast_df = blast_df_gene_dict[
                            genome_1_id + "_" + genome_2_id + ".blast.txt"
                        ]
                    else:
                        blast_df = blast_df_gene_dict[
                            genome_2_id + "_" + genome_1_id + ".blast.txt"
                        ]

                    # Case 1: Complete neighborhood OR aligned neighborhood with identical contig ends
                    if len(clust_1["genes"]) == len(clust_2["genes"]):
                        # Case 1.a: Assume neither are reversed
                        for gene_1, gene_2 in zip(
                            range(len(clust_1["genes"])), range(len(clust_2["genes"]))
                        ):
                            if not equivalent:
                                break
                            # Get locus tags
                            genome_1_qid = clust_1["genes"][gene_1]["uid"]
                            genome_2_qid = clust_2["genes"][gene_2]["uid"]

                            # Get gene names
                            genome_1_gene = clust_1["genes"][gene_1]["name"]
                            genome_2_gene = clust_2["genes"][gene_2]["name"]

                            filtered_df = blast_df.filter(
                                (
                                    (blast_df["query_id"] == genome_1_qid)
                                    & (blast_df["sub_id"] == genome_2_qid)
                                )
                                | (blast_df["query_id"] == genome_2_qid)
                                & (blast_df["sub_id"] == genome_1_qid)
                            )

                            # Check if match found (by locus tag or gene names)
                            if (
                                filtered_df.is_empty()
                                or filtered_df.filter(pl.col("PI") >= 70).is_empty()
                            ) and genome_1_gene != genome_2_gene:
                                # Case 1.b: Try reversing second neighborhood
                                for gene_1, gene_2 in zip(
                                    range(len(clust_1["genes"])),
                                    reversed(range(len(clust_2["genes"]))),
                                ):
                                    if not equivalent:
                                        break
                                    # Get locus tags
                                    genome_1_qid = clust_1["genes"][gene_1]["uid"]
                                    genome_2_qid = clust_2["genes"][gene_2]["uid"]

                                    # Get gene names
                                    genome_1_gene = clust_1["genes"][gene_1]["name"]
                                    genome_2_gene = clust_2["genes"][gene_2]["name"]

                                    filtered_df = blast_df.filter(
                                        (
                                            (blast_df["query_id"] == genome_1_qid)
                                            & (blast_df["sub_id"] == genome_2_qid)
                                        )
                                        | (
                                            (
                                                (blast_df["query_id"] == genome_2_qid)
                                                & (blast_df["sub_id"] == genome_1_qid)
                                            )
                                        )
                                    )

                                    if (
                                        filtered_df.is_empty()
                                        and genome_1_gene != genome_2_gene
                                    ):
                                        equivalent = False
                                        break
                                    elif (
                                        filtered_df.filter(
                                            pl.col("PI") >= 70
                                        ).is_empty()
                                        and genome_1_gene != genome_2_gene
                                    ):
                                        equivalent = False
                                        break

                    # Case 2: At least one neighborhood incomplete, ends don't align
                    else:
                        # Determine the smaller and larger neighborhoods
                        smaller_genome = (
                            clust_1
                            if len(clust_1["genes"]) < len(clust_2["genes"])
                            else clust_2
                        )
                        smaller_neighborhood = (
                            clust_1["genes"]
                            if len(clust_1["genes"]) < len(clust_2["genes"])
                            else clust_2["genes"]
                        )
                        larger_neighborhood = (
                            clust_2["genes"]
                            if len(clust_1["genes"]) < len(clust_2["genes"])
                            else clust_1["genes"]
                        )

                        # Determine the indices of the focal gene in the smaller and larger neighborhoods
                        with open(
                            output_path + "/neighborhood_indices.json", "r"
                        ) as infile:
                            if len(infile.readlines()) != 0:
                                infile.seek(0)
                                neighborhood_indices_json = json.load(infile)

                        for gene_name in neighborhood_indices_json:
                            if gene_name == extract_gene:
                                genome_indices = neighborhood_indices_json[gene_name][
                                    clust_1["name"]
                                ]
                                focal_gene_index_smaller = get_focal_index(
                                    genome_indices, num_neighbors
                                )

                                # Align the smaller neighborhood to the larger neighborhood using the focal genes
                                aligned_smaller = smaller_neighborhood[
                                    focal_gene_index_smaller : focal_gene_index_smaller
                                    + len(larger_neighborhood)
                                ]
                                aligned_smaller_reversed = smaller_neighborhood[
                                    focal_gene_index_smaller : focal_gene_index_smaller
                                    + len(larger_neighborhood)
                                ][::-1]
                                # If empty, not equivalent!
                                if not aligned_smaller and not aligned_smaller_reversed:
                                    equivalent = False
                                    break

                                # Iterate over the aligned genes and compare them in pairs
                                for gene_index, (gene, gene_reversed) in enumerate(
                                    zip(aligned_smaller, aligned_smaller_reversed)
                                ):
                                    if not equivalent:
                                        break
                                    # Get locus tags
                                    gene_1 = gene["uid"]
                                    gene_2 = larger_neighborhood[gene_index]["uid"]
                                    gene_1_reversed = gene_reversed["uid"]

                                    # Get gene names
                                    genome_1_gene = gene["name"]
                                    genome_2_gene = larger_neighborhood[gene_index][
                                        "name"
                                    ]
                                    genome_1_gene_rev = gene_reversed["name"]

                                    filtered_df = blast_df.filter(
                                        (
                                            (blast_df["query_id"] == gene_1)
                                            & (blast_df["sub_id"] == gene_2)
                                        )
                                        | (
                                            (blast_df["query_id"] == gene_2)
                                            & (blast_df["sub_id"] == gene_1)
                                        )
                                    )
                                    filtered_df_reversed = blast_df.filter(
                                        (
                                            (blast_df["query_id"] == gene_1_reversed)
                                            & (blast_df["sub_id"] == gene_2)
                                        )
                                        | (
                                            (blast_df["query_id"] == gene_2)
                                            & (blast_df["sub_id"] == gene_1_reversed)
                                        )
                                    )

                                    if (
                                        filtered_df.is_empty()
                                        and filtered_df_reversed.is_empty()
                                        and genome_1_gene != genome_2_gene
                                        and genome_1_gene_rev != genome_2_gene
                                    ):
                                        equivalent = False
                                        break
                                    elif (
                                        filtered_df.filter(
                                            pl.col("PI") >= 70
                                        ).is_empty()
                                        and filtered_df_reversed.filter(
                                            pl.col("PI") >= 70
                                        ).is_empty()
                                        and genome_1_gene != genome_2_gene
                                        and genome_1_gene_rev != genome_2_gene
                                    ):
                                        equivalent = False
                                        break

                    if equivalent:
                        represented_genomes.append(genome_2)
                        genome_genes_dict[genome_2] = True

        representative_genomes[genome] = represented_genomes

    # Ensure we don't represent a genome as a surrogate for itself
    delete_keys = []
    for genome, represented_genomes in representative_genomes.items():
        if genome_genes_dict[genome]:
            delete_keys.append(genome)

    for key in delete_keys:
        del representative_genomes[key]

    return representative_genomes


def write_filtered_genomes_textfile(representative_genomes, gene, output_path):
    """
    Writes textfile containing list of all surrogate genomes separated by a colon followed by the list of genomes they
    are representing (i.e. which have identical conserved neighborhoods).
    """
    out_path = output_path + "/JSON/surrogates/"
    check_output_path(out_path)

    with open(out_path + gene + "_surrogates.txt", "w+") as outfile:
        for representative in representative_genomes.keys():
            # Check if representative is unique or not
            if representative_genomes[representative]:
                outfile.write(representative + ": ")
                outfile.write(
                    str(representative_genomes[representative]).replace("'", "") + "\n"
                )


def filter_Loose_hits(RGI_dataframes_dict, PI=0.80):
    """
    Returns RGI dataframes with only Loose hits above or equal to a given threshold of Percent Identity.
    If user does not specify PI for Best_Identities as cutoff, the default value is 80%.
    """
    for genome_id, RGI_df in RGI_dataframes_dict.items():
        # Retrieve all Loose hits for the given genome
        try:
            filtered_df = RGI_df[
                ~((RGI_df["Cut_Off"] == "Loose") and (RGI_df["Best_Identities"] < PI))
            ]
            RGI_dataframes_dict[genome_id] = filtered_df
        except KeyError:
            print("No Loose hits found in the genome!")

    return RGI_dataframes_dict
