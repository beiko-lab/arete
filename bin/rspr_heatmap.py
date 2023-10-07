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
    return parser.parse_args(args)


def make_heatmap(results, output_path):
    print("Generating heatmap")
    data = (
        results.groupby(["tree_size", "exact_drSPR"]).size().reset_index(name="count")
    )
    data_pivot = data.pivot(
        index="exact_drSPR", columns="tree_size", values="count"
    ).fillna(0)
    sns.heatmap(
        data_pivot.loc[sorted(data_pivot.index, reverse=True)], annot=True, fmt=".0f"
    ).set(title="Number of trees")
    plt.xlabel("Tree size")
    plt.ylabel("Exact rSPR distance")
    plt.savefig(output_path)


def make_heatmap_from_csv(input_path, output_path):
    print("Generating heatmap from CSV")
    results = pd.read_csv(input_path)
    make_heatmap(results, output_path)


def main(args=None):
    args = parse_args(args)

    results = pd.read_csv(args.DF)
    make_heatmap(results, args.OUTPUT)


if __name__ == "__main__":
    sys.exit(main())
