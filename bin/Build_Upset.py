#!/usr/bin/env python

import argparse
import pandas as pd
from matplotlib import cm
from matplotlib import pyplot as plt
from upsetplot import UpSet, from_memberships
import warnings

# Suppress the UserWarning for tight_layout
warnings.simplefilter(action='ignore', category=FutureWarning)

def load_data(samplesheet_file, feature_profile_file):
    """
    Load samplesheet and feature profile data from TSV files.

    Parameters:
        samplesheet_file (str): Path to the samplesheet file.
        feature_profile_file (str): Path to the feature profile file.

    Returns:
        pd.DataFrame: Samplesheet dataframe.
        pd.DataFrame: Feature profile dataframe.
    """
    try:
        samplesheet_df = pd.read_csv(samplesheet_file)
        feature_df = pd.read_csv(feature_profile_file, sep='\t')
    except FileNotFoundError as e:
        print(f"Error: {e}")
        raise

    # Convert all columns in samplesheet_df to strings
    samplesheet_df = samplesheet_df.astype(str)

    return samplesheet_df, feature_df

def process_data(samplesheet_df, feature_df, column_name):
    """
    Process the data to merge samplesheet and feature dataframes and create presence/absence matrix.

    Parameters:
        samplesheet_df (pd.DataFrame): Samplesheet dataframe.
        feature_df (pd.DataFrame): Feature profile dataframe.
        column_name (str): Column name for the feature.

    Returns:
        pd.DataFrame: Transposed presence/absence dataframe with groups column.
    """
    merged_df = feature_df.merge(samplesheet_df, on='genome', how='left')

    if column_name not in merged_df.columns:
        print(f"Warning: Column '{column_name}' not found in samplesheet. Skipping.")
        return None

    if merged_df[column_name].nunique() <= 1:
        print(f"Warning: Column '{column_name}' contains only a single unique value. Skipping.")
        return None

    selected_columns = ['genome', column_name] + list(feature_df.columns[1:])
    feature_df = merged_df[selected_columns]
    presence = feature_df.groupby(column_name).any().astype(int)
    presence.drop('genome', axis=1, inplace=True, errors='ignore')
    presence_transposed = presence.transpose()
    presence_transposed.reset_index(inplace=True)
    presence_transposed.rename(columns={'index': 'feature'}, inplace=True)
    presence_transposed['category'] = presence_transposed['feature'].str.split('_', expand=True)[0]

    presence_transposed['groups'] = presence_transposed.apply(
        lambda row: ','.join([col for col in presence_transposed.columns[2:] if row[col] == 1]), axis=1)

    presence_transposed = presence_transposed.loc[:, ['feature', 'category', 'groups']]

    # print("Presence Transposed DataFrame with 'groups' column:\n", presence_transposed.head())  # Debug print to check the dataframe state

    return presence_transposed

def create_upset_plot(samplesheet_file, feature_profile_file, column_names, output_file_prefix):
    """
    Create an Upset plot from samplesheet and feature profile files.

    Parameters:
        samplesheet_file (str): Path to the samplesheet file.
        feature_profile_file (str): Path to the feature profile file.
        column_names (list of str): List of column names for features.
        output_file_prefix (str): Prefix for the output Upset plot files.
    """
    samplesheet_df, feature_df = load_data(samplesheet_file, feature_profile_file)

    for column_name in column_names:
        presence_transposed = process_data(samplesheet_df, feature_df, column_name)

        if presence_transposed is None:
            continue

        if 'groups' not in presence_transposed.columns:
            print(f"Error: 'groups' column was not created for {column_name}. Skipping.")
            continue

        # Ensure 'groups' column is processed correctly
        memberships = presence_transposed['groups'].str.split(",").apply(lambda x: frozenset(x))
        rows_by_group = from_memberships(memberships, data=presence_transposed.set_index('feature'))

        # print("Rows by group:\n", rows_by_group.head())  # Debug print to check rows_by_group structure

        # Create UpSet plot
        fig, ax = plt.subplots(figsize=(24, 12))  # Increased width for better layout
        upset = UpSet(rows_by_group, show_counts=True, intersection_plot_elements=0, min_subset_size=1)
        upset.add_stacked_bars(by="category", colors=cm.Pastel1, title="category", elements=10)
        upset.plot(fig=fig)

        # Remove the black frame by turning off the top and right spines
        ax.spines[['top','right','left','bottom']].set_visible(False)
        ax.tick_params(left=False)
        ax.tick_params(bottom=False)
        ax.tick_params(labelleft=False)
        ax.tick_params(labelbottom=False)

        # Adjust layout to move legend outside the figure
        plt.tight_layout(rect=[0.1, 0.1, 0.75, 1])  # Adjust rect to leave space on the right for legend
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

        output_file = f"{output_file_prefix}_{column_name}.png"
        # Overwrite the existing file if it exists
        plt.savefig(output_file)
        plt.close(fig)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create Upset plot from samplesheet and feature profile files.')
    parser.add_argument('-s', '--samplesheet_file', type=str, required=True, help='Path to the samplesheet file in TSV format')
    parser.add_argument('-f', '--feature_profile_file', type=str, required=True, help='Path to the feature profile file in TSV format')
    parser.add_argument('-c', '--columns', type=str, required=True, help='Comma-separated list of column names for features')
    parser.add_argument('-o', '--output_file_prefix', type=str, required=True, help='Prefix for the output Upset plot files (e.g., upset_plot)')
    args = parser.parse_args()

    column_names = args.columns.split(',')

    create_upset_plot(args.samplesheet_file, args.feature_profile_file, column_names, args.output_file_prefix)
