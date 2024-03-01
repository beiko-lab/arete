#!/usr/bin/env python

# Create heatmap of similarity thresholds based on PopPUNK distances

from matplotlib import pyplot as plt


from io import StringIO
import sys
import subprocess
import itertools as it
import pandas as pd
import seaborn as sns
import argparse
import sys


def parse_args(args=None):
    Description = "Subset genomes and create heatmap of similarity thresholds based on PopPUNK distances."
    Epilog = "Example usage: python filter_genomes.py <FILE_IN> <FILE_OUT> -c <CORE_DISTANCE> -a <ACC_DISTANCE>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input distances file.")
    parser.add_argument("FILE_OUT", help="Output plot.")
    parser.add_argument(
        "-c",
        "--core_threshold",
        dest="CORE_THRESH",
        type=float,
        help="Core genome distance threshold",
    )
    parser.add_argument(
        "-a",
        "--accessory_threshold",
        dest="ACC_THRESH",
        type=float,
        help="Accessory distance threshold",
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
        _, to_remove = filter_genomes(file_in, core_thresh, acc_thresh)
        n_to_remove = len(to_remove)
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
    run_cmd = [
        "awk",
        "-v",
        f"core_dist={core}",
        "-v",
        f"acc_dist={acc}",
        '{ if ($3 <= core_dist && $4 <= acc_dist) { if (!($1 in checked)) { print $1 "," $2; checked[$2] = 1 } } }',
        str(file_in),
    ]

    proc_out = subprocess.Popen(
        run_cmd,
        stdout=subprocess.PIPE,
    )

    b = StringIO(proc_out.communicate()[0].decode("utf-8"))

    genome_mapping = pd.read_csv(b, sep=",", names=["kept", "removed"])

    removed = set(genome_mapping["removed"].to_list())

    return genome_mapping, removed


def write_removed_genomes(file_in, core, acc):
    genome_mapping, removed = filter_genomes(file_in, core, acc)

    genome_mapping.to_csv("genome_mapping.csv", index=False)

    with open("removed_genomes.txt", "w") as f:
        f.write("\n".join(set(removed)))


def main(args=None):
    args = parse_args(args)
    make_heatmap(args.FILE_IN, args.FILE_OUT)
    write_removed_genomes(args.FILE_IN, args.CORE_THRESH, args.ACC_THRESH)


if __name__ == "__main__":
    sys.exit(main())
