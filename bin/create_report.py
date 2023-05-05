#!/usr/bin/env python

import os
import sys
import gzip
import argparse
from Bio import SeqIO
from pandas import read_table, merge, DataFrame, isna, NA
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
        "-r", "--vfdb_fasta", dest="VFDB_FASTA", help="VFDB reference FASTA."
    )
    parser.add_argument(
        "-m",
        "--mobsuite_out",
        dest="MOBSUITE",
        help="Mob Recon outputs.",
        nargs="?",
        const=None,
    )
    return parser.parse_args(args)


def summarize_alignment(path, db_name):
    df = read_table(path)

    summary = df.copy()[["genome_id", "qseqid", "sseqid", "pident", "qcover"]]

    summary["qcover"] = summary["qcover"].where(summary["qcover"] <= 1.0, 1.0)
    summary["qcover"] = round(summary["qcover"], 2) * 100

    summary = summary.rename(
        columns={
            "qseqid": "orf",
            "sseqid": db_name,
            "pident": f"{db_name}_identity",
            "qcover": f"{db_name}_qcover",
        }
    )

    return summary


def parse_vfdb_fasta(vfdb_fasta):
    with gzip.open(vfdb_fasta, "rt") as handle:
        record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    vfdb_df = DataFrame.from_dict(
        {k: v.description for (k, v) in record_dict.items()},
        orient="index",
        columns=["vfdb_desc"],
    )

    return vfdb_df


def create_vfdb_report(df, vfdb_df):
    w_vfdb = df.copy().merge(vfdb_df, left_on="vfdb", right_index=True, how="left")

    w_vfdb[["vfdb_short_id", "vfdb1", "vfdb2"]] = (
        w_vfdb["vfdb_desc"].str.extractall("\(([^()]\w+\/?\w+)\)").unstack()
    )

    w_vfdb["vfdb"] = w_vfdb[["vfdb", "vfdb1", "vfdb2"]].apply(
        lambda row: "/".join(row.values.astype(str)) if not isna(row.values[0]) else NA,
        axis=1,
    )
    w_vfdb = w_vfdb.drop(columns=["vfdb_desc", "vfdb1", "vfdb2"])

    return w_vfdb


def create_report(ann, diamond_outs, rgi, vfdb_fasta, mobsuite):
    # Summarize DIAMOND outs
    diamond_sums = [
        summarize_alignment(out, os.path.basename(out).strip(".txt").lower())
        for out in diamond_outs
    ]

    # Get VFDB df
    vfdb_df = parse_vfdb_fasta(vfdb_fasta)

    # RGI output
    rgi_df = read_table(rgi)
    rgi_sum = rgi_df[rgi_df["Best_Identities"] > 80]
    rgi_sum = rgi_sum[["Contig", "Best_Hit_ARO", "Cut_Off", "genome_id"]].rename(
        columns={"Contig": "orf", "Best_Hit_ARO": "AMR", "Cut_Off": "rgi_cutoff"}
    )
    rgi_sum["orf"] = rgi_sum["orf"].str.rsplit("_", n=1).str.get(0)

    diamond_sums.append(rgi_sum)

    # Bakta/Prokka output
    ann_tool = os.path.basename(ann).strip(".txt").lower()

    ann_df = read_table(ann)

    if ann_tool == "bakta":
        ann_sum = ann_df[
            ["genome_id", "#Sequence Id", "Start", "Stop", "Locus Tag"]
        ].rename(columns={"#Sequence Id": "contig_id", "Locus Tag": "orf"})
    else:
        ann_sum = ann_df[["genome_id", "locus_tag", "length_bp"]].rename(
            columns={"locus_tag": "orf"}
        )

    ann_sum = ann_sum[~ann_sum["orf"].isnull()]

    # Merge results
    orf_based_merged = reduce(
        lambda left, right: merge(left, right, on=["genome_id", "orf"], how="outer"),
        diamond_sums,
    )

    orf_ann = ann_sum.merge(orf_based_merged, on=["genome_id", "orf"], how="inner")

    w_vfdb = create_vfdb_report(orf_ann, vfdb_df)

    if mobsuite is not None and ann_tool == "bakta":
        # MobRecon output
        mobrecon = read_table(mobsuite)
        mobrecon_plasmids = mobrecon[mobrecon["molecule_type"] == "plasmid"]
        mobrecon_sum = mobrecon_plasmids[
            ["sample_id", "contig_id", "primary_cluster_id"]
        ].rename(columns={"sample_id": "genome_id", "primary_cluster_id": "plasmid"})
        mobrecon_sum["contig_id"] = mobrecon_sum["contig_id"].str.extract("(contig\d+)")
        mobrecon_sum["contig_id"] = mobrecon_sum["contig_id"].str.replace(
            r"(?<=g)0+", "_"
        )

        mobsuite_ann = ann_sum.merge(
            mobrecon_sum, on=["genome_id", "contig_id"], how="inner"
        )

        merged_full = mobsuite_ann.merge(
            w_vfdb, on=["genome_id", "orf", "contig_id", "Start", "Stop"], how="outer"
        )

        merged_full.to_csv(
            path_or_buf="annotation_report.tsv.gz", sep="\t", index=False
        )

    else:
        w_vfdb.to_csv(path_or_buf="annotation_report.tsv.gz", sep="\t", index=False)


def main(args=None):
    args = parse_args(args)
    create_report(args.ANN, args.DIAMOND_OUTS, args.RGI, args.VFDB_FASTA, args.MOBSUITE)


if __name__ == "__main__":
    sys.exit(main())
