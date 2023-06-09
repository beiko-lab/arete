#!/usr/bin/env python

"""
Utils for file manipulation and general data preprocessing used by multiple modules.
"""

import os
import glob
import string
import random

def strip_brackets(var):
    """
    Preprocesses list or string with brackets to be bracketless string (e.g. for easier manipulation of dataframes).
    """
    if isinstance(var, list):
        final = str(var)
        clean_str = final.strip('[]')
    else:
        clean_str = var.strip('[]')

    return clean_str


def get_full_filepaths(dir_path):

    filenames = os.listdir(dir_path)
    filepaths = [dir_path + '/' + filenames[i] for i in range(len(filenames))]

    return filepaths


def get_dir_filepaths(filepaths, file_ext):
    """
    Given a path to the folder where all genome AMR neighborhoods for a given AMR gene are stored, returns a list
    of all neighborhood fasta filepaths.
    """
    try:
        neighborhood_filepaths = glob.glob(os.path.join(filepaths, file_ext))
        return neighborhood_filepaths
    except FileNotFoundError:
        print("Error: no {} files were found in the specified directory.".format(file_ext))
        
def get_fasta_filename(filepath):
    """
    Gets filename from filepath (i.e. removes directories and file suffix).
    """
    filename = os.path.basename(filepath).split(".fasta")[0]
    return filename


def get_filename(filepath, extract=False):
    """
    Gets filename from filepath (i.e. removes directories and file suffix).
    """
    if extract:
        if '.fna' in filepath:
            filename = os.path.basename(filepath).split('.')[0]
        else:
            filename = os.path.basename(filepath).split('_')[0]
    else:
        filename = os.path.basename(filepath).split(".")[0]
    return filename


def get_dir_filenames(filepaths):
    """
    Given list of filepaths, returns their names (i.e. filename stripped of suffix).
    """
    filenames = []
    for file_path in filepaths:
        filenames.append(os.path.basename(file_path).split(".")[0])

    return filenames


def check_output_path(path):
    """
    Checks for presence of specified output path and creates the directory if it does not exist
    """
    if not os.path.exists(path):
        os.makedirs(path)


def make_output_file(output_path):
    """
    Creates an empty textfile in the specified output directory.
    """
    with open(output_path, 'w') as file:
        pass
    file.close()


def move_file(new_path, current_path):
    """
    Moves file from src to dest directory.
    """
    full_new_path = os.path.join(current_path, new_path)
    os.rename(current_path, full_new_path)


def remove_files(path_to_dir, file_ext):
    """
    Removes all files from a specified directory that do not have the given extension.
    """
    for file in os.listdir(path_to_dir):
        if not file.endswith(file_ext):
            os.remove(path_to_dir + '/' + file)


def generate_alphanumeric_string(length):
    """
    Generate a random alphanumeric string of fixed length
    """
    str = string.ascii_lowercase
    return ''.join(random.choice(str) for i in range(length))


def make_fasta_contig_dict(fasta_path, gene):
    """
    For a given gene and given FASTA neighborhood files, creates a dictionary where keys correspond to indices from
    0 to N, and values correspond to contig ids.
    """
    fasta_dict = dict()
    fasta_files = [fasta_file for fasta_file in os.listdir(fasta_path + '/' + gene)
                   if fasta_file.endswith('.fasta')]

    for fasta_file in fasta_files:
        with open(fasta_path + '/' + gene + '/' + fasta_file, 'r') as infile:
            data = infile.readlines()
            fasta_file_dict = dict()
            index = 0
            contig_lines = [line for line in data if line.startswith('>')]
            for line in contig_lines:
                contig_id = line.strip().replace('>', '')
                fasta_file_dict[index] = contig_id
                index += 1
            genome_id = fasta_file.split('.fasta')[0]
            fasta_dict[genome_id] = fasta_file_dict

    return fasta_dict