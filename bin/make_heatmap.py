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
    parser.add_argument(
        "-c",
        "--core_threshold",
        dest="CORE_THRESH",
        type=float,
        help="Core genome threshold",
    )
    parser.add_argument(
        "-a",
        "--accessory_threshold",
        dest="ACC_THRESH",
        type=float,
        help="Accessory threshold",
    )
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
        to_remove = filter_genomes(file_in, core_thresh, acc_thresh)
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
        title="Number of unique genomes by similarity threshold value"
    )
    plt.savefig(file_out)


def filter_genomes(file_in, core, acc):
    df = pd.read_table(file_in)

    checked = []
    # Iterate through the DataFrame
    for _, row in df.iterrows():
        query_genome = row["Query"]
        reference_genome = row["Reference"]
        core_value = row["Core"]
        accessory_value = row["Accessory"]

        if reference_genome in checked:
            continue
        else:
            if core_value <= core and accessory_value <= acc:
                checked.append(query_genome)

    return checked

def write_removed_genomes(file_in, core, acc):
    removed = filter_genomes(file_in, core, acc)

    with open("removed_genomes.txt", "w") as f:
        f.write("\n".join(set(removed)))


def main(args=None):
    args = parse_args(args)
    make_heatmap(args.FILE_IN, args.FILE_OUT)
    write_removed_genomes(args.FILE_IN, args.CORE_THRESH, args.ACC_THRESH)


if __name__ == "__main__":
    sys.exit(main())
