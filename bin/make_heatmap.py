#!/usr/bin/env python

# Create heatmap of similarity thresholds based on PopPUNK distances

from matplotlib import pyplot as plt
import itertools as it
import pandas as pd
import seaborn as sns
import argparse
import sys


def parse_args(args=None):
    Description = "Create heatmap of similarity thresholds based on PopPUNK distances."
    Epilog = "Example usage: python make_heatmap.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input distances file.")
    parser.add_argument("FILE_OUT", help="Output plot.")
    return parser.parse_args(args)


def make_heatmap(file_in, file_out):
    # Distance table
    df = pd.read_table(file_in)
    total_genomes = len(set(df["Query"].to_list() + df["Reference"].to_list()))
    # define the thresholds to investigate
    thresholds = [0.9, 0.95, 0.99, 0.999, 0.9999, 1.0]
    records = []

    for t in it.product(thresholds, thresholds):
        core_thresh = round(1 - t[0], 6)
        acc_thresh = round(1 - t[1], 6)
        # Which genomes can be removed
        to_remove = df[(df["Core"] <= core_thresh) & (df["Accessory"] <= acc_thresh)][
            "Query"
        ].to_list()
        n_to_remove = len(set(to_remove))
        n_genomes = total_genomes - n_to_remove

        records.append(
            {
                "core": round(t[0], 4),
                "accessory": round(t[1], 4),
                "n_genomes": n_genomes,
            }
        )

    # make dataframe from the info
    ndf = pd.DataFrame.from_records(records)
    ndf = ndf.pivot(index="core", columns="accessory", values="n_genomes")
    # plot
    sns.heatmap(ndf.loc[sorted(ndf.index, reverse=True)], annot=True, fmt=".0f").set(
        title="Number of unique genomes by threshold value"
    )
    plt.savefig(file_out)


def main(args=None):
    args = parse_args(args)
    make_heatmap(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
