#!/usr/bin/env python

"""
Given RGI files and their corresponding GBK files for complete bacterial genomes, as well as gene neighborhood clusters
generated using neighborhood_clustering.py, the utilities in this class are used for two main types of filtering:

a) Filtering of redundant AMR gene neighborhoods based on the clustering results for cleaner visualization downstream.
b) If user opts to investigate Loose hits using predictive module, filtering of low percent identity Loose hits.
"""

from utils import check_output_path

def check_sublist(neighborhood_1_genes, neighborhood_2_genes):
    length = len(neighborhood_2_genes)
    for i in range(len(neighborhood_1_genes) - length + 1):
        if all(neighborhood_2_genes[j] == neighborhood_1_genes[i+j] for j in range(length)):
            return True
    return False


def filter_neighborhoods(neighborhoods_dict, num_neighbors):
    """
    Given a dictionary containing genome identifiers as keys as neighborhood dataframes as values for an AMR gene, finds
    all identical neighborhood representations and keeps only one to render as a representative sample in the final gene
    order visualization.
    Outputs [AMR_gene]_representatives.txt file listing which genomes the representative neighborhood is also
    standing in for.
    """
    genome_genes_dict = {}
    for genome, neighborhood_df in neighborhoods_dict.items():
        genes = neighborhood_df['Gene_Name'].tolist()
        for g in range(len(genes)):
            gene = genes[g].strip('"')
            if gene.startswith('UID'):
                gene = 'UID'
            genes[g] = gene
        genome_genes_dict[genome] = [genes, False]

    # Store representative genomes in dict as: {'representative_genome': ['genomeA', 'genomeB', 'genomeC']}
    representative_genomes = {}
    candidate_rep_genomes = list(genome_genes_dict.keys())
    genomes = list(genome_genes_dict.keys())
    for genome in candidate_rep_genomes:
        if genome not in representative_genomes.keys():
            represented_genomes = []

            if genome_genes_dict[genome][1] == False:
                for genome_2 in genomes:

                    # If neighborhoods are identical
                    reversed_genome_2 = genome_genes_dict[genome_2][0][::-1]
                    is_contig_end_fwd = check_sublist(genome, genome_2)
                    is_contig_end_rev = check_sublist(genome, reversed_genome_2)

                    # Check if it's identical (forward or reversed)
                    if (genome != genome_2) and ((genome_genes_dict[genome][0] == genome_genes_dict[genome_2][0])
                                                 or (genome_genes_dict[genome][0] == reversed_genome_2)):
                        represented_genomes.append(genome_2)
                        genome_genes_dict[genome_2][1] = True
                    #elif (is_contig_end_fwd or is_contig_end_rev):
                    #    represented_genomes.append(genome_2)
                    #    genome_genes_dict[genome_2][1] = True

                representative_genomes[genome] = represented_genomes

    # Ensure we don't represent a genome as a surrogate for itself
    delete_keys_1 = []
    for genome, represented_genomes in representative_genomes.items():
        if len(represented_genomes) == 1 and represented_genomes[0] == genome:
            delete_keys_1.append(genome)


    delete_keys_2 = []
    for genome, represented_genomes in representative_genomes.items():
        for i in range(len(represented_genomes)-1):
            if represented_genomes[i] == genome:
                del represented_genomes[i]

    for key in delete_keys_1:
        del representative_genomes[key]
    #for key_2 in delete_keys_2:
    #    del representative_genomes[]

    return representative_genomes


def write_filtered_genomes_textfile(representative_genomes, AMR_gene, output_path):
    """
    Writes textfile containing list of all surrogate genomes separated by a colon followed by the list of genomes they
    are representing (i.e. which have identical conserved neighborhoods).
    """
    out_path = output_path + '/JSON/surrogates/'
    check_output_path(out_path)

    with open(out_path + AMR_gene + '_surrogates.txt', 'w+') as outfile:
        for representative in representative_genomes.keys():
            # Check if representative is unique or not
            if representative_genomes[representative] == []:
                outfile.write(representative + '\n')
            else:
                outfile.write(representative + ': ')
                outfile.write(str(representative_genomes[representative]).replace("'", "") + '\n')


def filter_Loose_hits(RGI_dataframes_dict, PI=0.80):
    """
    Returns RGI dataframes with only Loose hits above or equal to a given threshold of Percent Identity.
    If user does not specify PI for Best_Identities as cutoff, the default value is 80%.
    """
    for genome_id, RGI_df in RGI_dataframes_dict.items():
        # Retrieve all Loose hits for the given genome
        try:
            filtered_df = RGI_df[~((RGI_df['Cut_Off'] == 'Loose') and (RGI_df['Best_Identities'] < PI))]
            RGI_dataframes_dict[genome_id] = filtered_df
        except KeyError:
            print("No Loose hits found in the genome!")

    return RGI_dataframes_dict
