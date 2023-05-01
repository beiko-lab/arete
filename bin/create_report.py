#!/usr/bin/env python

import os
import sys
import argparse
from pandas import read_table, merge
from functools import reduce


def parse_args(args=None):
    Description = "Filter alignment results."
    Epilog = "Example usage: python filter_alignment.py <ANN> <DIAMOND_OUTS> <RGI> <MOBSUITE>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-a",
        "--annotation_out",
        dest="ANN",
        help="Annotation report (Bakta or Prokka).",
    )
    parser.add_argument(
        "-d",
        "--diamond_outs",
        dest="DIAMOND_OUTS",
        help="DIAMOND alignment outputs.",
        nargs="*",
    )
    parser.add_argument("-r", "--rgi_out", dest="RGI", help="RGI output.")
    parser.add_argument(
        "-m", "--mobsuite_out", dest="MOBSUITE", help="Mob Recon outputs."
    )
    return parser.parse_args(args)


def summarize_alignment(path, db_name):
    df = read_table(path)

    summary = df[["genome_id", "qseqid", "sseqid", "pident"]]

    summary = summary.rename(
        columns={"qseqid": "orf", "sseqid": db_name, "pident": f"{db_name}_identity"}
    )

    return summary


def create_report(ann, diamond_outs, rgi, mobsuite):
    # Summarize DIAMOND outs
    diamond_sums = [
        summarize_alignment(out, os.path.basename(out).strip(".txt").lower())
        for out in diamond_outs
    ]

    # RGI output
    rgi_df = read_table(rgi)
    rgi_sum = rgi_df[["Contig", "Best_Hit_ARO", "Cut_Off", "genome_id"]].rename(
        columns={"Contig": "orf", "Best_Hit_ARO": "AMR", "Cut_Off": "rgi_cutoff"}
    )
    rgi_sum["orf"] = rgi_sum["orf"].str.rsplit("_", n=1).str.get(0)

    orf_based_anns = diamond_sums.append(rgi_sum)

    # Bakta/Prokka output
    ann_df = read_table(ann)
    ann_sum = ann_df[
        ["genome_id", "#Sequence Id", "Start", "Stop", "Locus Tag"]
    ].rename(columns={"#Sequence Id": "contig_id", "Locus Tag": "orf"})
    ann_sum = ann_sum[~ann_sum["orf"].isnull()]

    # MobRecon output
    mobrecon = read_table(mobsuite)
    mobrecon_plasmids = mobrecon[mobrecon["molecule_type"] == "plasmid"]
    mobrecon_sum = mobrecon_plasmids[
        ["sample_id", "contig_id", "primary_cluster_id"]
    ].rename(columns={"sample_id": "genome_id", "primary_cluster_id": "plasmid"})
    mobrecon_sum["contig_id"] = mobrecon_sum["contig_id"].str.extract("(contig\d+)")
    mobrecon_sum["contig_id"] = mobrecon_sum["contig_id"].str.replace(r"(?<=g)0+", "_")

    # Merge results
    orf_based_merged = reduce(
        lambda left, right: merge(left, right, on=["genome_id", "orf"], how="outer"),
        orf_based_anns,
    )

    mobsuite_ann = ann_sum.merge(
        mobrecon_sum, on=["genome_id", "contig_id"], how="inner"
    )
    orf_ann = ann_sum.merge(orf_based_merged, on=["genome_id", "orf"], how="inner")

    merged_full = mobsuite_ann.merge(
        orf_ann, on=["genome_id", "orf", "contig_id", "Start", "Stop"], how="outer"
    )
    merged_full.to_csv(path_or_buf="annotation_report.tsv", sep="\t", index=False)


def main(args=None):
    args = parse_args(args)
    create_report(args.ANN, args.DIAMOND_OUTS, args.RGI, args.MOBSUITE)


if __name__ == "__main__":
    sys.exit(main())
