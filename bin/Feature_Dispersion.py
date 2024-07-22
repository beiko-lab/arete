#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from ete3 import Tree

def round_to_sig_figs(value, sig_figs):
    if value == 0:
        return 0
    return round(value, sig_figs - int(np.floor(np.log10(abs(value)))) - 1)

def calculate_phylogenetic_diversity(tree):
    sum_branch_lengths = 0.0
    for node in tree.traverse():
        sum_branch_lengths += node.dist
    return sum_branch_lengths

def calculate_feature_counts(feature_file):
    feature_df = pd.read_csv(feature_file, sep='\t', index_col=0)
    feature_df.reset_index(inplace=True)
    feature_df.rename(columns={'index': 'genome_id'}, inplace=True)

    # Initialize lists to store data for the new dataframe
    features = []
    total_counts = []
    genomes_list = []

    # Iterate over each feature column
    for feature_col in feature_df.columns[1:]:  # Exclude 'genome_id' column
        # Collect the genomes where this feature is present into a comma-separated list
        genomes_with_feature = feature_df[feature_df[feature_col] == 1]['genome_id'].tolist()
        genomes_str = ','.join(genomes_with_feature)

        # Calculate the total count of this feature across all genomes
        total_count_feature = feature_df[feature_col].sum()

        # Append data to lists
        features.append(feature_col)
        total_counts.append(total_count_feature)
        genomes_list.append(genomes_str)

    # Create a new dataframe from the collected data
    features_df = pd.DataFrame({
        'feature': features,
        'total_count': total_counts,
        'genomes_list': genomes_list
    })

    return features_df

def verify_genome_ids(tree, feature_file, samplesheet_file=None):
    feature_df = pd.read_csv(feature_file, sep='\t', index_col=0)
    feature_genomes = set(feature_df.index)

    if samplesheet_file:
        samplesheet_df = pd.read_csv(samplesheet_file, usecols=[0], header=0, names=['genome_id'], skiprows=1)
        samplesheet_genomes = set(samplesheet_df['genome_id'])
    else:
        samplesheet_genomes = set()

    tree_genomes = set(tree.get_leaf_names())

    missing_in_tree = feature_genomes.union(samplesheet_genomes) - tree_genomes

    if missing_in_tree:
        print(f"Error: The following genome IDs are missing in the phylogenetic tree: {', '.join(missing_in_tree)}")
        exit(1)

def generate_heatmap(output_df, output_heatmap):
    bins = np.arange(0, 1.1, 0.1)
    max_genome_count = output_df['Genome Count'].max()
    genome_bins = np.linspace(0, max_genome_count, 11)  # Generate 11 edges to create 10 bins
    genome_bins_labels = [f'{int(genome_bins[i]) + 1}-{int(genome_bins[i + 1])}' for i in range(len(genome_bins) - 1)]

    heatmap_data = pd.DataFrame(0, index=genome_bins_labels, columns=bins)

    for index, row in output_df.iterrows():
        pd_ratio = row['PD Ratio']
        genome_count = row['Genome Count']
        bin_idx = np.digitize(pd_ratio, bins) - 1
        genome_bin_idx = np.digitize(genome_count, genome_bins) - 1

        # Ensure the indices are within the valid range
        bin_idx = min(bin_idx, len(bins) - 1)
        genome_bin_idx = min(genome_bin_idx, len(genome_bins_labels) - 1)

        heatmap_data.iloc[genome_bin_idx, bin_idx] += 1

    plt.figure(figsize=(12, 8))
    sns.heatmap(np.log1p(heatmap_data), cmap="Reds", cbar_kws={'label': 'Number of Features (log scale)'}, annot=heatmap_data, fmt='g', linewidths=.5)
    plt.xlabel('PD Ratio Bins')
    plt.ylabel('Genome Count Bins')
    plt.title('Heatmap of Features by PD Ratio and Genome Count')
    plt.xticks(ticks=np.arange(0.5, len(bins), 1), labels=np.round(bins, 1))
    plt.yticks(ticks=np.arange(0.5, len(genome_bins_labels), 1), labels=genome_bins_labels)
    plt.gca().invert_yaxis()
    plt.savefig(output_heatmap)
    plt.close()





def main(tree_file, feature_file, output_base, samplesheet_file=None, samplesheet_columns=None):
    ref_tree = Tree(tree_file)

    # Verify genome IDs
    verify_genome_ids(ref_tree, feature_file, samplesheet_file)

    # Calculate phylogenetic diversity
    total_diversity = calculate_phylogenetic_diversity(ref_tree)

    # Calculate feature counts
    feature_distr = calculate_feature_counts(feature_file)

    # Read samplesheet if provided
    if samplesheet_file:
        samplesheet_df = pd.read_csv(samplesheet_file, header=0)
        available_columns = set(samplesheet_df.columns)
        if 'genome_id' not in available_columns:
            print("Error: 'genome_id' column is missing in the samplesheet.")
            exit(1)

        samplesheet_data = {}
        for column in samplesheet_columns:
            if column in available_columns:
                samplesheet_data[column] = samplesheet_df[['genome_id', column]].set_index('genome_id')[column].to_dict()
            else:
                print(f"Warning: Column '{column}' not found in the samplesheet. Skipping.")
    else:
        samplesheet_data = {}

    # Create an empty DataFrame to store the output
    output_columns = [
        'Feature Name',
        'Total PD', 'Projected PD', 'PD Ratio',
        'Genome Count','PD Ratio / Genome Count'
    ]

    # Add columns for each requested samplesheet column
    for column in samplesheet_columns:
        output_columns.append(f'{column} Distinct Values')
        output_columns.append(f'PD Ratio / {column} Values')

    output_df = pd.DataFrame(columns=output_columns)

    # Iterate over each feature in feature_distr
    for index, row in feature_distr.iterrows():
        feature_name = row['feature']
        genomes_list = row['genomes_list'].split(',')
        genome_count = row['total_count']

        # Only proceed if the feature is present in more than one genome
        if len(genomes_list) > 1:
            # Generate a list of genomes to keep (those that have the feature)
            genomes_to_keep = [genome for genome in genomes_list if genome in ref_tree]

            # Create a copy of the original tree with only the relevant genomes
            projected_tree = ref_tree.copy()
            projected_tree.prune(genomes_to_keep)

            # Calculate phylogenetic diversity of the projected tree
            projected_diversity = calculate_phylogenetic_diversity(projected_tree)

            # Calculate the ratio of projected diversity to total diversity
            ratio_diversity = projected_diversity / total_diversity

            # Calculate the ratio of projected phylogenetic diversity to total count of genomes
            genome_ratio_phylogenetic_diversity = ratio_diversity / genome_count

            # Prepare the row for output
            output_row = [
                feature_name,
                round_to_sig_figs(total_diversity, 4),
                round_to_sig_figs(projected_diversity, 4),
                round_to_sig_figs(ratio_diversity, 4),
                round_to_sig_figs(genome_count, 4),
                round_to_sig_figs(genome_ratio_phylogenetic_diversity, 4)
            ]

            # Add values for each requested samplesheet column
            for column in samplesheet_columns:
                if column in samplesheet_data:
                    # Identify distinct values for the column
                    distinct_values = set(samplesheet_data[column].get(genome, None) for genome in genomes_list if genome in samplesheet_data[column])
                    distinct_values.discard(None)
                    V = len(distinct_values)
                    PD_ratio_per_V = ratio_diversity / V if V > 0 else 0
                    output_row.extend([V, round_to_sig_figs(PD_ratio_per_V, 4)])
                else:
                    output_row.extend([None, None])

            # Add the row to the output DataFrame
            output_df.loc[index] = output_row

    output_df_sorted = output_df.sort_values(by='PD Ratio / Genome Count')

    # Save the output dataframe to a TSV file
    output_df_sorted.to_csv(f"{output_base}.tsv", sep='\t', index=False)

    # Generate the heatmap
    generate_heatmap(output_df_sorted, f"{output_base}.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate feature statistics based on phylogenetic tree and genus information.')
    parser.add_argument('--tree_file', type=str, required=True, help='Path to the Newick tree file')
    parser.add_argument('--feature_file', type=str, required=True, help='Path to the feature presence/absence file')
    parser.add_argument('--output_base', type=str, required=True, help='Base name for the output files (without extension)')
    parser.add_argument('--samplesheet_file', type=str, help='Path to the file mapping genome IDs to other properties')
    parser.add_argument('--samplesheet_columns', type=str, help='Comma-separated list of columns to process from the samplesheet')
    args = parser.parse_args()

    if args.samplesheet_columns:
        samplesheet_columns = args.samplesheet_columns.split(',')
    else:
        samplesheet_columns = []

    main(args.tree_file, args.feature_file, args.output_base, args.samplesheet_file, samplesheet_columns)
