#!/usr/bin/env python

# Run rspr

import sys
import os
import re
from pathlib import Path
import argparse
import subprocess
from ete3 import Tree
import pandas as pd
from collections import defaultdict
from matplotlib import pyplot as plt
import seaborn as sns


#####################################################################
### FUNCTION PARSE_ARGS
### Parse the command-line arguments.
### args: the arguments
### RETURN a list that contains values for all the defined arguments.
#####################################################################


def parse_args(args=None):
    Description = "Run rspr"
    Epilog = "Example usage: rspr.py CORE_TREE GENE_TREES"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("-c", "--core", dest="CORE_TREE", help="Core tree")
    parser.add_argument(
        "-a", "--acc", dest="GENE_TREES", help="Gene tree samplesheet path"
    )
    parser.add_argument(
        "-ann",
        "--annotation",
        dest="ANNOTATION",
        help="Annotation table from BAKTA/PROKKA",
        nargs="?",
    )
    parser.add_argument(
        "-o", "--output", dest="OUTPUT_DIR", default="approx", help="Gene tree list"
    )
    parser.add_argument(
        "-mrd",
        "--min_rspr_distance",
        dest="MIN_RSPR_DISTANCE",
        type=int,
        default=10,
        help="Minimum rSPR distance used to define processing groups.",
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


def read_tree(input_path):
    with open(input_path, "r") as f:
        tree_string = f.read()
        formatted = re.sub(r";[^:]+:", ":", tree_string)
        return Tree(formatted)


#####################################################################
### FUNCTION ROOT_TREE
### Root the unrooted input trees
### input_path: path of the input tree
### output_path: path of the output rooted trees
### RETURN rooted tree and tree size
#####################################################################


def root_tree(input_path, output_path):
    tre = read_tree(input_path)
    midpoint = tre.get_midpoint_outgroup()
    tre.set_outgroup(midpoint)
    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))
    tre.write(outfile=output_path)
    return tre.write(), len(tre.get_leaves())


#####################################################################
### FUNCTION ROOT_TREE
### Root all the unrooted input trees in directory
### core_tree: path of the core tree
### gene_trees: path of the csv file containing all the gene tree paths
### output_dir: output directory path
### results: dataframe of the results to store tree size
### merge_pair: boolean to check whether to merge coer tree and gene tree in a single file
### RETURN path of the rooted gene trees directory
#####################################################################


def root_trees(core_tree, gene_trees_path, output_dir, results, merge_pair=False):
    print("Rooting trees")
    #'''
    reference_tree = core_tree

    Path(output_dir, "rooted_reference_tree").mkdir(exist_ok=True, parents=True)
    Path(output_dir, "rooted_gene_trees").mkdir(exist_ok=True)

    rooted_reference_tree = os.path.join(
        output_dir, "rooted_reference_tree/core_gene_alignment.tre"
    )
    refer_content, refer_tree_size = root_tree(reference_tree, rooted_reference_tree)

    df_gene_trees = pd.read_csv(gene_trees_path)
    rooted_gene_trees_path = os.path.join(output_dir, "rooted_gene_trees")
    for filename in df_gene_trees["path"]:
        basename = Path(filename).name
        rooted_gene_tree_path = os.path.join(rooted_gene_trees_path, basename)
        gene_content, gene_tree_size = root_tree(filename, rooted_gene_tree_path)
        results.loc[basename, "tree_size"] = gene_tree_size
        if merge_pair:
            with open(rooted_gene_tree_path, "w") as f2:
                f2.write(refer_content + "\n" + gene_content)
    #'''
    return rooted_gene_trees_path


#####################################################################
### FUNCTION EXTRACT_APPROX_DISTANCE
### Extract approx rspr distance from the rspr output
### text: rspr output text
### RETURN approx rspr distance
#####################################################################


def extract_approx_distance(text):
    for line in text.splitlines():
        if "approx drSPR=" in line:
            distance = line.split("approx drSPR=")[1].strip()
            return distance
    return "0"


#####################################################################
### FUNCTION APPROX_RSPR
### Run approx rspr algorithm of set of input tree pairs
### rooted_gene_trees_path: path of the rooted gene trees directory
### results: dataframe to store the approx rspr results
### min_branch_len: minimum branch length
### max_support_threshold: maximum branching support threshold
#####################################################################


def approx_rspr(
    rooted_gene_trees_path, results, min_branch_len=0, max_support_threshold=0.7
):
    print("Calculating approx distance")
    rspr_path = [
        "rspr",
        "-approx",
        "-multifurcating",
        "-length " + str(min_branch_len),
        "-support " + str(max_support_threshold),
    ]
    for filename in os.listdir(rooted_gene_trees_path):
        gene_tree_path = os.path.join(rooted_gene_trees_path, filename)
        with open(gene_tree_path, "r") as infile:
            result = subprocess.run(
                rspr_path, stdin=infile, capture_output=True, text=True
            )
            dist = extract_approx_distance(result.stdout)
            results.loc[filename, "approx_drSPR"] = dist


#####################################################################
### FUNCTION GENERATE_HEATMAP
### Generate heatmap figure from frequency table
### freq_table: frequency table of tree size and approx rspr distance
### output_path: output path for storing the heatmap
#####################################################################


def generate_heatmap(freq_table, output_path):
    print("Generating heatmap")
    sns.heatmap(
        freq_table, annot=True, fmt=".0f"
    ).set(title="Number of trees")
    plt.xlabel("Tree size")
    plt.ylabel("Approx rSPR distance")
    plt.savefig(output_path)
    plt.clf()

#####################################################################
### FUNCTION MAKE_HEATMAP
### Generate heatmap of tree size and approx rspr distance
### results: dataframe of the approx rspr distance and tree size
### output_path: output path for storing the heatmap
#####################################################################

def make_heatmap(results, output_path):
    print("Generating heatmap")
    data = (
        results.groupby(["tree_size", "approx_drSPR"]).size().reset_index(name="count")
    )
    data_pivot = data.pivot(
        index="approx_drSPR", columns="tree_size", values="count"
    ).fillna(0)
    generate_heatmap(data_pivot.loc[sorted(data_pivot.index, reverse=True)], output_path)


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
### Generate heatmap of tree size and approx rspr distance groups
### results: dataframe of the approx rspr distance and tree size groups
### output_path: output path for storing the heatmap
#####################################################################

def make_group_heatmap(results, output_path):
    print("Generating group heatmap")
    data = pd.crosstab(results["approx_drSPR"], results["tree_size"])

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

    plt.figure(figsize=(12, 12))
    generate_heatmap(aggregated_row_df, output_path)


#####################################################################
### FUNCTION GENERATE_GROUP_SIZES
### Generate groups sizes based on number of tree pairs avaialble
### target_sum: total number of trees available
### RETURN groups of trees
#####################################################################

def generate_group_sizes(target_sum, max_groups=500):
    degree = 1
    current_sum = 0
    group_sizes = []
    while current_sum < target_sum:
        k = 0
        current_sum = 0
        group_sizes = []
        while k < max_groups:
            value = int(2 ** (k /(max_groups / degree)))
            if value > 1e20:
                break
            group_sizes.append(value)
            current_sum += value
            k += 1
        degree += 1
    return group_sizes


#####################################################################
### FUNCTION MAKE_GROUPS
### Generate groups of tree pairs based on the approx rspr distnace
### results: the input results dataframe
### RETURN groups of trees
#####################################################################

def make_groups_v1(results, min_limit=10):
    print("Generating groups")
    min_group = results[results["approx_drSPR"] <= min_limit]["file_name"].tolist()
    groups = defaultdict()
    first_group = "group_0"
    groups[first_group] = min_group

    rem_results = results[results["approx_drSPR"] > min_limit].sort_values(
        by="approx_drSPR", ascending=False
    )
    rem_length = len(rem_results)
    group_sizes = generate_group_sizes(rem_length)
    cur_index, grp_idx = 0, 0
    while cur_index < rem_length:
        cur_group_names = rem_results.iloc[cur_index : cur_index + group_sizes[grp_idx]]["file_name"].tolist()
        groups[f"group_{grp_idx+1}"] = cur_group_names
        cur_index += group_sizes[grp_idx]
        grp_idx += 1
    return groups


#####################################################################
### FUNCTION MAKE_GROUPS
### Generate groups of tree pairs based on the approx rspr distnace
### min_limit: minimum approx rspr distance for the first group
### RETURN groups of trees
#####################################################################

def make_groups(results, min_limit=10):
    print("Generating groups")
    min_group = results[results["approx_drSPR"] <= min_limit]["file_name"].tolist()
    groups = defaultdict()
    first_group = "group_0"
    groups[first_group] = min_group

    rem_results = results[results["approx_drSPR"] > min_limit].sort_values(
        by="approx_drSPR", ascending=False
    )
    rem_length = len(rem_results)
    cur_index, grp_size, cur_appx_dist, grp_name = 0, 0, -1, 1
    while cur_index < rem_length:
        if cur_appx_dist != rem_results.iloc[cur_index]["approx_drSPR"]:
            cur_appx_dist = rem_results.iloc[cur_index]["approx_drSPR"]
            grp_size += 1
        cur_group_names = rem_results.iloc[cur_index : cur_index + grp_size][
            "file_name"
        ].tolist()
        groups[f"group_{grp_name}"] = cur_group_names
        cur_index += grp_size
        grp_name += 1
    return groups


def make_groups_from_csv(input_df, min_limit):
    print("Generating groups from CSV")
    groups = make_groups_v1(input_df, min_limit)
    tidy_data = [
        (key, val)
        for key, value in groups.items()
        for val in (value if isinstance(value, list) else [value])
    ]
    df_w_groups = pd.DataFrame(tidy_data, columns=["group_name", "file_name"])
    merged = input_df.merge(df_w_groups, on="file_name")

    return merged


def make_heatmap_from_csv(input_path, output_path):
    print("Generating heatmap from CSV")
    results = pd.read_csv(input_path)
    make_heatmap(results, output_path)


def join_annotation_data(df, annotation_data):
    ann_df = pd.read_table(annotation_data, dtype={"genome_id": "str"})
    ann_df.columns = map(str.lower, ann_df.columns)
    ann_df.columns = ann_df.columns.str.replace(" ", "_")
    ann_subset = ann_df[["gene", "product"]]

    df["tree_name"] = [f.split(".")[0] for f in df["file_name"]]

    merged = df.merge(ann_subset, how="left", left_on="tree_name", right_on="gene")

    if merged["gene"].isnull().all():
        ann_subset = ann_df[["locus_tag", "gene", "product"]]
        merged = df.merge(
            ann_subset, how="left", left_on="tree_name", right_on="locus_tag"
        )

    return merged.fillna(value="NULL").drop("tree_name", axis=1).drop_duplicates()


def main(args=None):
    args = parse_args(args)

    # Approx RSPR
    #'''
    results = pd.DataFrame(columns=["file_name", "tree_size", "approx_drSPR"])
    results.set_index("file_name", inplace=True)
    rooted_paths = root_trees(
        args.CORE_TREE, args.GENE_TREES, args.OUTPUT_DIR, results, True
    )
    approx_rspr(
        rooted_paths,
        results,
        args.MIN_BRANCH_LENGTH,
        args.MAX_SUPPORT_THRESHOLD,
    )
    results["approx_drSPR"] = pd.to_numeric(results["approx_drSPR"])
    fig_path = os.path.join(args.OUTPUT_DIR, "output.png")
    make_heatmap(results, fig_path)
    group_fig_path = os.path.join(args.OUTPUT_DIR, "group_output.png")
    make_group_heatmap(results, group_fig_path)

    results.reset_index("file_name", inplace=True)
    if args.ANNOTATION:
        results = join_annotation_data(results, args.ANNOTATION)
    res_path = os.path.join(args.OUTPUT_DIR, "output.tsv")
    df_with_groups = make_groups_from_csv(results, args.MIN_RSPR_DISTANCE)
    df_with_groups.to_csv(res_path, sep="\t", index=False)
    groups = df_with_groups["group_name"].unique()
    grouped_dfs = [
        df_with_groups[df_with_groups["group_name"] == group] for group in groups
    ]
    for df in grouped_dfs:
        group = df["group_name"].unique()[0]
        df.to_csv(f"rspr_{group}.csv", index=False)

    #'''

    # From CSV
    """
    phylo_dir = os.path.join(args.INPUT_DIR_PATH, "phylogenomics")
    res_path = os.path.join(phylo_dir, 'output.csv')
    fig_path = os.path.join(phylo_dir, 'output.png')
    make_heatmap_from_csv(res_path, fig_path)
    """


if __name__ == "__main__":
    sys.exit(main())
