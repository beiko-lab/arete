#!/usr/bin/env python

# Run rspr

import sys
import os
import argparse
import subprocess
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict
from ete3 import Tree, TreeStyle
import re

#####################################################################
### FUNCTION PARSE_ARGS
### Parse the command-line arguments.
### args: the arguments
### RETURN a list that contains values for all the defined arguments.
#####################################################################


def parse_args(args=None):
    Description = "Run rspr"
    Epilog = "Example usage: rspr.py INPUT_DIR_PATH"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-df", "--dataframe", dest="SUBSET_DF", help="CSV table subset by group"
    )
    parser.add_argument(
        "-mar",
        "--max_approx_rspr",
        dest="MAX_APPROX_RSPR_DIST",
        type=int,
        default=-1,
        help="Maximum approximate rSPR distance",
    )
    parser.add_argument(
        "-mbl",
        "--min_branch_length",
        dest="MIN_BRANCH_LENGTH",
        type=int,
        default=0,
        help="Minimum branch length",
    )
    parser.add_argument(
        "-mst",
        "--max_support_threshold",
        dest="MAX_SUPPORT_THRESHOLD",
        type=float,
        default=0.7,
        help="Maximum support threshold",
    )
    return parser.parse_args(args)


#####################################################################
### FUNCTION EXTRACT_EXACT_DISTANCE
### Extract exact rspr distance from the rspr output
### text: rspr output text
### RETURN exact rspr distance
#####################################################################

def extract_exact_distance(text):
    for line in text.splitlines():
        if "total exact drSPR=" in line:
            distance = line.split("total exact drSPR=")[1].strip()
            return distance
    return "0"


#####################################################################
### FUNCTION GET_CLUSTER_BRANCH_VAL
### Calculate cluster network branch value
### lst_child: list of child of node
### dict_clstr_map: cluster map
#####################################################################

def get_cluster_branch_val(lst_child, dict_clstr_map):
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
### FUNCTION UPDATE_BRANCH_LENGTHS_TO_CLUSTER_VAL
### Update branch length of the cluster trees to cluster value
### node: current node
### dict_clstr_map: cluster map
#####################################################################

def update_branch_lengths_to_cluster_val(node, dict_clstr_map, total_trees, leaf_mapping):
    lst_node = list()
    if not node:
        return lst_node

    if node.is_leaf():
        lst_node.append(leaf_mapping[node.name])
        node.dist = (get_cluster_branch_val(lst_node, dict_clstr_map) / total_trees) * 100
    else:
        for child in node.children:
            lst_child_node = update_branch_lengths_to_cluster_val(child, dict_clstr_map, total_trees, leaf_mapping)
            lst_node.extend(lst_child_node)
        node.dist = (get_cluster_branch_val(lst_node, dict_clstr_map) / total_trees) * 100
    return lst_node


#####################################################################
### FUNCTION GENERATE_CLUSTER_NETWORK
### Generate cluster network from list of clusters
### lst_tree_clusters: list of clusters for all gene trees
### refer_tree: reference tree
#####################################################################

def generate_cluster_network(lst_tree_clusters, refer_tree):
    print("Generating cluster network")
      
    leaf_mapping = {leaf: i for i, leaf in enumerate(refer_tree.get_leaf_names())}
    dict_clstr_map = defaultdict(int)
    for tree in lst_tree_clusters:
        for cluster in tree:
            dict_clstr_map[str(sorted(cluster))] += 1

    total_trees = len(lst_tree_clusters)
    update_branch_lengths_to_cluster_val(refer_tree.get_tree_root(), dict_clstr_map, total_trees, leaf_mapping)


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


#####################################################################
### FUNCTION FPT_RSPR
### Run exact rspr algorithm of set of input tree pairs
### results_df: dataframe to input tree pairs and store the exact rspr results
### min_branch_len: minimum branch length
### max_support_threshold: maximum branching support threshold
#####################################################################

def fpt_rspr(results_df, min_branch_len=0, max_support_threshold=0.7, gather_cluster_info=True):
    print("Calculating exact distance")
    current_dir = os.path.dirname(os.path.abspath(__file__))
    exe_path = os.path.join(current_dir, 'rspr.exe')
    rspr_path = [
        exe_path,
        "-multifurcating",
        "-show_clusters",
        "-length " + str(min_branch_len),
        "-support " + str(max_support_threshold),
    ]

    lst_tree_clusters = []
    trees_path = os.path.join("rooted_gene_trees")

    # Run this groups in parallel
    for filename in results_df.index:
        gene_tree_path = os.path.join(trees_path, filename)
        full_output = ""
        with open(gene_tree_path, "r") as infile:
            if not gather_cluster_info:
            # Clustering information not required
                result = subprocess.run(
                    rspr_path, stdin=infile, capture_output=True, text=True
                )
                full_output = result.stdout
            else:
            # Clustering information required
                process = subprocess.Popen(
                    rspr_path, stdin=infile, stdout=subprocess.PIPE, text=True, universal_newlines=True
                )

                clusters = []
                clustering_start = False
                output_lines = []
                while True:
                    line = process.stdout.readline()
                    if not line:
                        break
                    if "Clusters start" in line:
                        clustering_start = True
                        continue
                    elif "Clusters end" in line:
                        clustering_start = False
                    
                    if clustering_start:
                        updated_line = line.replace('(', '').replace(')', '').replace('\n', '')
                        cluster_nodes = updated_line.split(',')
                        cluster_nodes = [int(node) for node in cluster_nodes if "X" not in node]
                        clusters.append(cluster_nodes)
                    
                    output_lines.append(line)
                lst_tree_clusters.append(clusters)
                process.wait()
                full_output = ''.join(output_lines)

            dist = extract_exact_distance(full_output)
            results_df.loc[filename, "exact_drSPR"] = dist
    return lst_tree_clusters


def read_tree(input_path):
    with open(input_path, "r") as f:
        tree_string = f.read()
        formatted = re.sub(r";[^:]+:", ":", tree_string)
        return Tree(formatted)


def main(args=None):
    args = parse_args(args)

    # Exact RSPR
    csv_path = os.path.join(args.SUBSET_DF)
    results = pd.read_csv(csv_path)
    if args.MAX_APPROX_RSPR_DIST >= 0:
        results = results[(results["approx_drSPR"] <= args.MAX_APPROX_RSPR_DIST)]

    results.set_index("file_name", inplace=True)
    lst_tree_clusters = fpt_rspr(results, args.MIN_BRANCH_LENGTH, args.MAX_SUPPORT_THRESHOLD, True)

    refer_tree_path = os.path.join("rooted_reference_tree/core_gene_alignment.tre")
    refer_tree = read_tree(refer_tree_path)
    generate_cluster_network(lst_tree_clusters, refer_tree)

    cluster_tree_path = os.path.join("cluster_tree.png")
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    refer_tree.render(cluster_tree_path, tree_style=ts)

    group = results["group_name"].unique()[0]
    res_path = f"exact_output_{group}.tsv"
    results.fillna(value="NULL").to_csv(res_path, sep="\t", index=True)

    # From CSV
    """
    phylo_dir = os.path.join(args.INPUT_DIR_PATH, "phylogenomics")
    res_path = os.path.join(phylo_dir, 'exact_output.csv')
    fig_path = os.path.join(phylo_dir, 'exact_output.png')
    make_heatmap_from_csv(res_path, fig_path)
    """


if __name__ == "__main__":
    sys.exit(main())
