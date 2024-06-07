#!/usr/bin/env python

import re
import sys
import argparse
from pandas import read_table, read_csv, merge, DataFrame
from pathlib import Path


def parse_args(args=None):
    Description = "Create recombination workflow input data"
    Epilog = "Example usage: python organize_recomb_data.py <QUAST_REPORT> <POPPUNK_CLUSTERS> <ASSEMBLIES> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("QUAST_REPORT", help="QUAST QC report.")
    parser.add_argument(
        "POPPUNK_CLUSTERS", help="Genomes and the clusters they belong to."
    )
    parser.add_argument("ASSEMBLIES", help="Assembly samplesheet.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def get_cluster_samplesheet(df, cluster):
    out_path = Path(f"cluster_{cluster}.txt").resolve()
    subset = df.loc[df["Cluster"] == cluster, ["id", "path"]]
    subset["path"] = [Path(path).name for path in subset["path"].to_list()]

    subset.to_csv(out_path, sep="\t", index=False, header=False)

    return str(out_path)


def remove_after_second_underscore(val):
    parts = val.split('_')
    if len(parts) > 2:
        return '_'.join(parts[:2])
    return val

def create_recomb_input(quast_report, poppunk_clusters, assembly_samplesheet, file_out):
    # Parsing datasets
    quast = read_table(quast_report).loc[:, ["Assembly", "N50"]]
    quast["Assembly"] = (
            quast["Assembly"]
            .apply(remove_after_second_underscore)
            )

    poppunk = read_csv(poppunk_clusters)
    poppunk["Taxon"] = poppunk["Taxon"].str.replace("\.(.*)|_T1|$", "", regex=True)

    assemblies = read_csv(assembly_samplesheet, names=["id", "assembly_path"])
    assemblies["assembly_path"] = [
        Path(path).name for path in assemblies["assembly_path"].to_list()
    ]
    assemblies["id"] = assemblies["id"].str.replace("\.(.*)|_T1|$", "", regex=True)

    # Merging datasets
    merged = poppunk.merge(quast, right_on="Assembly", left_on="Taxon").loc[
        :, ["Cluster", "Assembly", "N50"]
    ]
    df_highest_n50 = (
        merged.sort_values("N50", ascending=False)
        .groupby("Cluster")
        .first()
        .reset_index()
    )
    df_highest_n50
    df_highest_n50.rename(columns={"Assembly": "Reference"}, inplace=True)
    df_merged = merge(
        merged, df_highest_n50[["Cluster", "Reference"]], on="Cluster", how="left"
    )
    df_with_assemblies = (
        merge(assemblies, df_merged, left_on="id", right_on="Assembly")
        .merge(
            assemblies,
            right_on="id",
            left_on="Reference",
        )
        .loc[
            :,
            ["Cluster", "Assembly", "assembly_path_x", "Reference", "assembly_path_y"],
        ]
        .rename(
            columns={
                "Assembly": "id",
                "assembly_path_x": "path",
                "assembly_path_y": "reference_path",
            }
        )
    )

    # Getting cluster samplesheets
    clusters = df_with_assemblies["Cluster"].unique()
    cluster_sheets = [
        get_cluster_samplesheet(df_with_assemblies, cluster) for cluster in clusters
    ]
    cluster_ids = [
        int(re.findall(r"cluster_(\d+)", cluster_path)[0])
        for cluster_path in cluster_sheets
    ]
    cluster_df = DataFrame({"Cluster": cluster_ids, "samplesheet": cluster_sheets})

    # Writing final dataframe
    final_df = cluster_df.merge(
        df_with_assemblies[["Cluster", "Reference", "reference_path"]], on="Cluster"
    ).drop_duplicates()

    final_df.to_csv(file_out, index=False)


def main(args=None):
    args = parse_args(args)
    create_recomb_input(
        args.QUAST_REPORT, args.POPPUNK_CLUSTERS, args.ASSEMBLIES, args.FILE_OUT
    )


if __name__ == "__main__":
    sys.exit(main())
