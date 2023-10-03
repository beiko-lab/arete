#!/usr/bin/env python

import sys
import argparse
import pandas as pd


def parse_args(args=None):
    Description = "Setup rspr groups"
    Epilog = "Example usage: setup_rspr_groups.py RSPR_DF"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("-df", "--dataframe", dest="DF", help="CSV table from Approx")
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    # Exact RSPR
    results = pd.read_csv(args.DF)
    groups = results["group_name"].unique()
    grouped_dfs = [results[results["group_name"] == group] for group in groups]
    for df in grouped_dfs:
        group = df["group_name"].unique()[0]
        df.to_csv(f"rspr_{group}.csv", index=False)
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
