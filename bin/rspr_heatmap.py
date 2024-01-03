#!/usr/bin/env python

import sys
import argparse
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
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

if __name__ == "__main__":
    sys.exit(main())
