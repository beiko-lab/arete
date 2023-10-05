#!/usr/bin/env python

# Run rspr

import sys
import os
import argparse
import subprocess
import pandas as pd
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
        type=int,
        default=0.7,
        help="Maximum support threshold",
    )
    return parser.parse_args(args)


# def make_heatmap(results, output_path):
#     print("Generating heatmap")
#     data = (
#         results.groupby(["tree_size", "exact_drSPR"]).size().reset_index(name="count")
#     )
#     data_pivot = data.pivot(
#         index="exact_drSPR", columns="tree_size", values="count"
#     ).fillna(0)
#     sns.heatmap(
#         data_pivot.loc[sorted(data_pivot.index, reverse=True)], annot=True, fmt=".0f"
#     ).set(title="Number of trees")
#     plt.xlabel("Tree size")
#     plt.ylabel("Exact rSPR distance")
#     plt.savefig(output_path)


# def make_heatmap_from_csv(input_path, output_path):
#     print("Generating heatmap from CSV")
#     results = pd.read_csv(input_path)
#     make_heatmap(results, output_path)


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
### FUNCTION FPT_RSPR
### Run exact rspr algorithm of set of input tree pairs
### results_df: dataframe to input tree pairs and store the exact rspr results
### min_branch_len: minimum branch length
### max_support_threshold: maximum branching support threshold
#####################################################################


def fpt_rspr(results_df, min_branch_len=0, max_support_threshold=0.7):
    print("Calculating exact distance")
    rspr_path = [
        "rspr",
        "-multifurcating",
        "-length " + str(min_branch_len),
        "-support " + str(max_support_threshold),
    ]
    trees_path = os.path.join("rooted_gene_trees")
    # Run this groups in parallel
    for filename in results_df.index:
        gene_tree_path = os.path.join(trees_path, filename)
        with open(gene_tree_path, "r") as infile:
            result = subprocess.run(
                rspr_path, stdin=infile, capture_output=True, text=True
            )
            dist = extract_exact_distance(result.stdout)
            results_df.loc[filename, "exact_drSPR"] = dist


def main(args=None):
    args = parse_args(args)

    # Exact RSPR
    results = pd.read_csv(args.SUBSET_DF)
    if args.MAX_APPROX_RSPR_DIST >= 0:
        results = results[(results["approx_drSPR"] <= args.MAX_APPROX_RSPR_DIST)]

    results.set_index("file_name", inplace=True)
    fpt_rspr(results, args.MIN_BRANCH_LENGTH, args.MAX_SUPPORT_THRESHOLD)

    group = results["group_name"].unique()[0]
    res_path = f"exact_output_{group}.csv"
    results.to_csv(res_path, index=True)
    # fig_path = os.path.join(phylo_dir, "exact_output.png")
    # make_heatmap(results, fig_path)

    # From CSV
    """
    phylo_dir = os.path.join(args.INPUT_DIR_PATH, "phylogenomics")
    res_path = os.path.join(phylo_dir, 'exact_output.csv')
    fig_path = os.path.join(phylo_dir, 'exact_output.png')
    make_heatmap_from_csv(res_path, fig_path)
    """


if __name__ == "__main__":
    sys.exit(main())
