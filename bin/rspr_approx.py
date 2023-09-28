#!/usr/bin/env python

# Run rspr

import sys
import os
import argparse
import subprocess
from ete3 import Tree
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

def parse_args(args=None):
    Description = "Run rspr"
    Epilog = "Example usage: rspr.py INPUT_DIR_PATH"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("INPUT_DIR_PATH", help="Input directory path")
    parser.add_argument("-mbl", "--min_branch_length", dest="MIN_BRANCH_LENGTH", type=int, default=0, help="Minimum branch length")
    parser.add_argument("-mst", "--max_support_threshold", dest="MAX_SUPPORT_THRESHOLD", type=int, default=0, help="Maximum support threshold")
    return parser.parse_args(args)

def root_tree(input_path, output_path):
    if os.path.isfile(input_path):
        tre = Tree(input_path)
        midpoint = tre.get_midpoint_outgroup()
        tre.set_outgroup(midpoint)
        if not os.path.exists(os.path.dirname(output_path)):
            os.makedirs(os.path.dirname(output_path))
        tre.write(outfile=output_path)
        return tre.write(), len(tre.get_leaves())

def root_trees(input_dir_path, results, merge_pair=False):
    print("Rooting trees")
    phylo_dir = os.path.join(input_dir_path, "phylogenomics")
    #'''
    reference_tree = os.path.join(phylo_dir, "reference_tree\core_gene_alignment.tre")
    rooted_reference_tree = os.path.join(phylo_dir, "rooted_reference_tree\core_gene_alignment.tre")
    refer_content, refer_tree_size = root_tree(reference_tree, rooted_reference_tree)

    gene_trees_path = os.path.join(phylo_dir, "gene_trees")
    rooted_gene_trees_path = os.path.join(phylo_dir, "rooted_gene_trees")
    for filename in os.listdir(gene_trees_path):
        gene_tree_path = os.path.join(gene_trees_path, filename)
        rooted_gene_tree_path = os.path.join(rooted_gene_trees_path, filename)
        gene_content, gene_tree_size = root_tree(gene_tree_path, rooted_gene_tree_path)
        results.loc[filename, "tree_size"] = gene_tree_size
        if merge_pair:
            with open(rooted_gene_tree_path, 'w') as f2:
                f2.write(refer_content + "\n" + gene_content)
    #'''
    return phylo_dir

def extract_approx_distance(text):
    for line in text.splitlines():
        if "approx drSPR=" in line:
            distance = line.split("approx drSPR=")[1].strip()
            return distance
    return "0"

def approx_rspr(input_path, results, min_branch_len=0, max_support_threshold=0):
    print("Calculating approx distance")
    current_dir = os.path.dirname(os.path.abspath(__file__))
    exe_path = os.path.join(current_dir, 'rspr.exe')
    rspr_path = [exe_path, "-approx", "-multifurcating", "-length " + str(min_branch_len), "-support " + str(max_support_threshold)]
    trees_path = os.path.join(input_path, "rooted_gene_trees")
    for filename in os.listdir(trees_path):
        gene_tree_path = os.path.join(trees_path, filename)
        with open(gene_tree_path, 'r') as infile:
            result = subprocess.run(rspr_path, stdin=infile, capture_output=True, text=True)
            dist = extract_approx_distance(result.stdout)
            results.loc[filename, "approx_drSPR"] = dist

def make_heatmap(results, output_path):
    print("Generating heatmap")
    data = results.groupby(["tree_size", "approx_drSPR"]).size().reset_index(name="count")
    data_pivot = data.pivot(index="approx_drSPR", columns="tree_size", values="count").fillna(0)
    sns.heatmap(data_pivot.loc[sorted(data_pivot.index, reverse=True)], annot=True, fmt=".0f").set(
        title="Number of trees"
    )
    plt.xlabel("Tree size")
    plt.ylabel("Approx rSPR distance")
    plt.savefig(output_path)

def make_heatmap_from_csv(input_path, output_path):
    print("Generating heatmap from CSV")
    results = pd.read_csv(input_path)
    make_heatmap(results, output_path)

def main(args=None):
    args = parse_args(args)

    # Approx RSPR
    #'''
    results = pd.DataFrame(columns=["file_name", "tree_size", "approx_drSPR"])
    results.set_index('file_name', inplace=True)
    phylo_dir = root_trees(args.INPUT_DIR_PATH, results, True)
    approx_rspr(phylo_dir, results, args.MIN_BRANCH_LENGTH, args.MAX_SUPPORT_THRESHOLD)

    res_path = os.path.join(phylo_dir, 'output.csv')
    results.to_csv(res_path, index=True)
    fig_path = os.path.join(phylo_dir, 'output.png')
    make_heatmap(results, fig_path)
    #'''

    # From CSV
    '''
    phylo_dir = os.path.join(args.INPUT_DIR_PATH, "phylogenomics")
    res_path = os.path.join(phylo_dir, 'output.csv')
    fig_path = os.path.join(phylo_dir, 'output.png')
    make_heatmap_from_csv(res_path, fig_path)
    '''


if __name__ == "__main__":
    sys.exit(main())