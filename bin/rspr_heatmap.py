#!/usr/bin/env python

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


def parse_args(args=None):
    Description = "Run rspr heatmap"
    Epilog = "Example usage: rspr_heatmap.py DF OUTPUT"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("-df", "--dataframe", dest="DF", help="CSV table")
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
    return parser.parse_args(args)


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
    plt.ylabel("Exact rSPR distance")
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
        results.groupby(["tree_size", "exact_drSPR"]).size().reset_index(name="count")
    )
    data_pivot = data.pivot(
        index="exact_drSPR", columns="tree_size", values="count"
    ).fillna(0)
    generate_heatmap(data_pivot.loc[sorted(data_pivot.index, reverse=True)], output_path)


#####################################################################
### FUNCTION GET_GROUP_SIZE
### Get preferred group size for generating heatmap
### all_tree_sizes: list of all values
### max_size: maximum number of groups
#####################################################################

def get_group_size(all_values, max_groups=15):
    group_size = 1
    if len(all_values) <= max_groups:
        return group_size

    multipliers = [2, 2.5, 2]
    cur_mul_idx = 0
    while not (all_values/group_size) <= max_groups:
        group_size *= multipliers[cur_mul_idx]
        cur_mul_idx = (cur_mul_idx + 1) % len(multipliers)
    return int(group_size)


#####################################################################
### FUNCTION MAKE_GROUP_HEATMAP
### Generate heatmap of tree size and exact rspr distance groups
### results: dataframe of the exact rspr distance and tree size groups
### output_path: output path for storing the heatmap
#####################################################################

def make_group_heatmap(results, output_path):
    print("Generating group heatmap")
    data = pd.crosstab(results["exact_drSPR"], results["tree_size"])

    all_tree_sizes = data.columns.astype('int32')
    tree_group_size = get_group_size(all_tree_sizes)
    aggregated_df = pd.DataFrame()
    if tree_group_size > 1:
        for i in range(1, max(all_tree_sizes), tree_group_size):
            group_columns = [col for col in all_tree_sizes if i <= int(col) <= i + tree_group_size - 1]
            group_sum = data[group_columns].sum(axis=1)
            aggregated_df[f'{i}-{i+tree_group_size-1}'] = group_sum
    else:
        aggregated_df = data

    all_distances = aggregated_df.index.astype('int32')
    distance_group_size = get_group_size(all_distances)
    aggregated_row_df = pd.DataFrame(columns=aggregated_df.columns)
    if distance_group_size > 1:
        for i in range(1, max(all_distances), distance_group_size):
            group_rows = [row for row in all_distances if i <= int(row) <= i + distance_group_size - 1]
            group_sum = aggregated_df.loc[group_rows].sum(axis=0)
            aggregated_row_df.loc[f'{i}-{i+distance_group_size-1}'] = group_sum
    else:
        aggregated_row_df = aggregated_df

    plt.figure(figsize=(14, 14))
    generate_heatmap(aggregated_row_df, output_path)


def make_heatmap_from_csv(input_path, output_path):
    print("Generating heatmap from CSV")
    results = pd.read_table(input_path)
    make_heatmap(results, output_path)


def main(args=None):
    args = parse_args(args)

    results = pd.read_table(args.DF)
    make_heatmap(results, args.OUTPUT)
    make_group_heatmap(results, args.GROUP_OUTPUT)

if __name__ == "__main__":
    sys.exit(main())
