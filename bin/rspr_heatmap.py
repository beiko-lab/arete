#!/usr/bin/env python

import sys
import argparse
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
import os 
import json
from Bio import Phylo
import re
import io
from collections import defaultdict

def parse_args(args=None):
    Description = "Run rspr heatmap"
    Epilog = "Example usage: rspr_heatmap.py DF OUTPUT"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("-df", "--dataframe", dest="DF", help="CSV table")
    parser.add_argument(
        "-cf",
        "--cluster_file",
        dest="CF",
        help="Cluster data file"
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="OUTPUT",
        help="Output file name",
    )
    parser.add_argument(
        "-go",
        "--group_output",
        dest="GROUP_OUTPUT",
        help="Group output file name",
    )
    parser.add_argument(
        "-co",
        "--cluster_output",
        dest="CLUSTER_OUTPUT",
        help="Cluster probability output file name",
    )
    parser.add_argument(
        "-cfo",
        "--cluster_file_output",
        dest="CLUSTER_FILE_OUTPUT",
        help="Cluster probability output newick tree file name",
    )
    parser.add_argument(
        "-mnher",
        "--min_heatmap_exact_rspr",
        dest="MIN_HEATMAP_RSPR_DISTANCE",
        type=int,
        default=0,
        help="Minimum exact rSPR distance used to generate heatmap.",
    )
    parser.add_argument(
        "-mxher",
        "--max_heatmap_exact_rspr",
        dest="MAX_HEATMAP_RSPR_DISTANCE",
        type=int,
        default=-1,
        help="Maximum exact rSPR distance used to generate heatmap.",
    )
    return parser.parse_args(args)


#####################################################################
### FUNCTION GENERATE_HEATMAP
### Generate heatmap figure from frequency table
### freq_table: frequency table of tree size and approx rspr distance
### output_path: output path for storing the heatmap
#####################################################################

def generate_heatmap(freq_table, output_path, log_scale=False):
    plt.figure(figsize=(12, 12))
    ax = sns.heatmap(
        freq_table, annot=True, fmt=".0f", norm=LogNorm() if log_scale else None
    )
    ax.invert_yaxis()
    plt.title("Number of trees")
    plt.xlabel("Tree size")
    plt.ylabel("Exact rSPR distance")
    plt.savefig(output_path)
    plt.clf()

#####################################################################
### FUNCTION MAKE_HEATMAP
### Generate heatmap of tree size and approx rspr distance
### results: dataframe of the approx rspr distance and tree size
### output_path: output path for storing the heatmap
#####################################################################

def make_heatmap(results, output_path, min_distance, max_distance):
    print("Generating heatmap")

    # create sub dataframe
    sub_results = results[(results["exact_drSPR"] >= min_distance)]
    if max_distance >= 0:
        sub_results = sub_results[(sub_results["exact_drSPR"] <= max_distance)]

    data = (
        results.groupby(["tree_size", "exact_drSPR"]).size().reset_index(name="count")
    )
    data_pivot = data.pivot(
        index="exact_drSPR", columns="tree_size", values="count"
    ).fillna(0)
    generate_heatmap(data_pivot.loc[sorted(data_pivot.index, reverse=True)], output_path)


def make_heatmap_from_tsv(input_path, output_path, min_distance, max_distance):
    print("Generating heatmap from CSV")
    results = pd.read_table(input_path)
    make_heatmap(results, output_path, min_distance, max_distance)

#####################################################################
### FUNCTION GET_GROUP_SIZE
### Get preferred group size for generating heatmap
### all_tree_sizes: list of all values
### max_size: maximum number of groups
#####################################################################

def get_heatmap_group_size(all_values, max_groups=15):
    group_size = 1
    if len(all_values) <= max_groups:
        return group_size

    multipliers = [2, 2.5, 2]
    cur_mul_idx = 0
    max_val = max(all_values)
    while not (max_val/group_size) <= max_groups:
        group_size *= multipliers[cur_mul_idx]
        cur_mul_idx = (cur_mul_idx + 1) % len(multipliers)
    return int(group_size)


#####################################################################
### FUNCTION MAKE_GROUP_HEATMAP
### Generate heatmap of tree size and exact rspr distance groups
### results: dataframe of the exact rspr distance and tree size groups
### output_path: output path for storing the heatmap
#####################################################################

def make_group_heatmap(results, output_path, min_distance, max_distance):
    print("Generating group heatmap")

    # create sub dataframe
    sub_results = results[(results["exact_drSPR"] >= min_distance)]
    if max_distance >= 0:
        sub_results = sub_results[(sub_results["exact_drSPR"] <= max_distance)]

    data = pd.crosstab(sub_results["exact_drSPR"], sub_results["tree_size"])

    all_tree_sizes = data.columns.astype('int32')
    tree_group_size = get_heatmap_group_size(all_tree_sizes)
    aggregated_df = pd.DataFrame()
    if tree_group_size > 1:
        for i in range(1, max(all_tree_sizes) + 1, tree_group_size):
            group_columns = [col for col in all_tree_sizes if i <= int(col) <= i + tree_group_size - 1]
            group_sum = data[group_columns].sum(axis=1)
            group_start = i if i > min(all_tree_sizes) else min(all_tree_sizes)
            group_end = (i+tree_group_size-1) if (i+tree_group_size-1) < max(all_tree_sizes) else max(all_tree_sizes)
            aggregated_df[f'{group_start}-{group_end}'] = group_sum
    else:
        aggregated_df = data

    all_distances = aggregated_df.index.astype('int32')
    distance_group_size = get_heatmap_group_size(all_distances)
    aggregated_row_df = pd.DataFrame(columns=aggregated_df.columns)
    if distance_group_size > 1:
        for i in range(0, max(all_distances) + 1, distance_group_size):
            group_rows = [row for row in all_distances if i <= int(row) <= i + distance_group_size - 1]
            group_sum = aggregated_df.loc[group_rows].sum(axis=0)
            group_start = i if i > min(all_distances) else min(all_distances)
            group_end = (i+distance_group_size-1) if (i+distance_group_size-1) < max(all_distances) else max(all_distances)
            aggregated_row_df.loc[f'{group_start}-{group_end}'] = group_sum
    else:
        aggregated_row_df = aggregated_df
    generate_heatmap(aggregated_row_df, output_path, True)


#region cluster network

#####################################################################
### FUNCTION GET_CLUSTER_PROBABILITY_VAL
### Calculate cluster probability value for node
### lst_child: list of child of node
### dict_clstr_map: cluster map
#####################################################################

def get_cluster_probability_val(lst_child, dict_clstr_map):
    lst_keys = []
    total_count = 0
    for key in dict_clstr_map:
        lst_key_child = [int(item) for item in key.strip('[]').split(', ')]
        if set(lst_key_child).issubset(set(lst_child)):
            total_count += dict_clstr_map[key]
            lst_keys.append(key)
    for key in lst_keys:
        del dict_clstr_map[key]
    return total_count


#####################################################################
### FUNCTION UPDATE_CLUSTER_PROBABILITY
### Update cluster probability value and assign it to node
### node: current node
### dict_clstr_map: cluster map
#####################################################################

def update_cluster_probability(node, dict_clstr_map, total_trees, leaf_mapping):
    lst_node = list()
    if not node:
        return lst_node

    if node.is_terminal():
        lst_node.append(leaf_mapping[node.name])
        node.cluster_probability = (get_cluster_probability_val(lst_node, dict_clstr_map) / total_trees) * 100
    else:
        for child in node.clades:
            lst_child_node = update_cluster_probability(child, dict_clstr_map, total_trees, leaf_mapping)
            lst_node.extend(lst_child_node)
        node.cluster_probability = (get_cluster_probability_val(lst_node, dict_clstr_map) / total_trees) * 100
    return lst_node


#####################################################################
### FUNCTION GENERATE_CLUSTER_NETWORK
### Generate cluster network from list of clusters
### lst_tree_clusters: list of clusters for all gene trees
### refer_tree: reference tree
#####################################################################

def generate_cluster_network(lst_tree_clusters, refer_tree):
    print("Generating cluster network")
    if not refer_tree:
        return
      
    lst_leaves = [leave.name for leave in refer_tree.get_terminals()]
    leaf_mapping = {leaf: i for i, leaf in enumerate(lst_leaves)}
    dict_clstr_map = defaultdict(int)
    for tree in lst_tree_clusters:
        for cluster in tree:
            if len(cluster) > 0:
                dict_clstr_map[str(sorted(cluster))] += 1

    total_trees = len(lst_tree_clusters)
    update_cluster_probability(refer_tree.root, dict_clstr_map, total_trees, leaf_mapping)


#####################################################################
### FUNCTION GENERATE_CLUSTER_HEATMAP
### Generate cluster heatmap
### lst_tree_clusters: list of clusters for all gene trees
### cluster_heatmap_path: putput path of heatmap
#####################################################################

def generate_cluster_heatmap(lst_tree_clusters, cluster_heatmap_path):
    print("Generating cluster heatmap")
    leaves = set(element for tree in lst_tree_clusters for cluster in tree for element in cluster)
    leaves = sorted(list(leaves))
    precentage_matrix = [[0 for _ in range(len(leaves))] for _ in range(len(leaves))]

    total_trees = len(lst_tree_clusters)
    for i, element1 in enumerate(leaves):
        for j, element2 in enumerate(leaves):
            if i < j:
                total_element_trees = 0
                element1_present = False
                element2_present = False
                count = 0
                for tree in lst_tree_clusters:
                    for cluster in tree:
                        if element1 in cluster:
                            element1_present = True
                        if element2 in cluster:
                            element2_present = True
                        if element1 in cluster and element2 in cluster:
                            count += 1
                            break
                    if element1_present and element2_present:
                        total_element_trees+=1
                percentage = (count / total_element_trees) * 100
                precentage_matrix[i][j] = percentage
                precentage_matrix[j][i] = percentage
            elif i == j:
                precentage_matrix[i][j] = 100

    plt.figure()
    sns.heatmap(precentage_matrix, annot=True, fmt=".0f").set(title="Cluster wise leaves distribution")
    plt.xlabel("Leaves")
    plt.ylabel("Leaves")
    plt.savefig(cluster_heatmap_path)

    
def read_tree(input_path):
    with open(input_path, "r") as f:
        tree_string = f.read()
        formatted = re.sub(r";[^:]+:", ":", tree_string)
        return Phylo.read(io.StringIO(formatted), "newick")

def write_tree(output_path, data):
    with open(output_path, "w") as f:
        f.write(data) 

def get_fig_size(refer_tree):
    max_fig_size = 100
    num_leaves = refer_tree.count_terminals()
    fig_size = num_leaves / 2
    if fig_size > max_fig_size:
        fig_size = max_fig_size
    return fig_size


#endregion
    
def main(args=None):
    args = parse_args(args)

    results = pd.read_table(args.DF)

    # Generate standard heatmap
    results["exact_drSPR"] = pd.to_numeric(results["exact_drSPR"])
    make_heatmap(
        results, 
        args.OUTPUT,
        args.MIN_HEATMAP_RSPR_DISTANCE,
        args.MAX_HEATMAP_RSPR_DISTANCE
    )

    # Generate group heatmap
    make_group_heatmap(
        results,
        args.GROUP_OUTPUT,
        args.MIN_HEATMAP_RSPR_DISTANCE,
        args.MAX_HEATMAP_RSPR_DISTANCE
    )

    # Generate cluster network
    lst_tree_clusters = []
    cluster_path = args.CF
    with open(cluster_path, "r") as file:
        for str_clstr in file:
            lst_tree_clusters.append(json.loads(str_clstr))

    cluster_tree_path = args.CLUSTER_OUTPUT
    cluster_file_path = args.CLUSTER_FILE_OUTPUT
    refer_tree_path = os.path.join("rooted_reference_tree/core_gene_alignment.tre")
    refer_tree = read_tree(refer_tree_path)
    if refer_tree:
        generate_cluster_network(lst_tree_clusters, refer_tree)
        write_tree(cluster_file_path, refer_tree)

        plt.rcParams['font.size'] = '12'
        fig_size = get_fig_size(refer_tree)
        fig, ax = plt.subplots(figsize=(fig_size, fig_size))
        Phylo.draw(refer_tree, axes=ax, do_show=False, branch_labels=lambda c: round(c.cluster_probability, 2))

        plt.xlabel('Cluster probability (%)')
        plt.title("Cluster network")
        plt.savefig(cluster_tree_path, format="png")

if __name__ == "__main__":
    sys.exit(main())
