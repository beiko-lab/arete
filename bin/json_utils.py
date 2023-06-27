#!/usr/bin/env python
"""
Methods for generating JSON representations of extracted neighborhoods that can be used with clustermap.js.
"""

import json
import itertools
import pandas as pd
from pathlib import Path
from utils import (
    check_output_path,
    generate_alphanumeric_string,
    make_fasta_contig_dict,
)


def swap_neighborhood_orientation(df):
    """
    Reverses neighborhood representation.
    Used to ensure all AMR genes have same orientation in downstream gene order visualizations.
    """
    df["Gene_Strand"] = df["Gene_Strand"].map({-1: +1, +1: -1}, na_action=None)

    return df


def reverse_start_end(df, n_start, n_stop):
    """
    Assuming other aspects of gene neighborhood df have been reversed (e.g. orientation, order, downstream/upstream),
    modifies each neighborhood gene's start and stop index to reflect these changes.
    """
    for gene in df.index:
        start = df["Gene_Start"][gene]
        end = df["Gene_End"][gene]

        # Subtract end coord
        adjusted_start = start - n_stop
        adjusted_end = end - n_stop

        # Flip values and swap
        new_start = abs(adjusted_end)
        new_end = abs(adjusted_start)

        # Update indices
        df.loc[gene, ["Gene_Start"]] = [new_start + n_start]
        df.loc[gene, ["Gene_End"]] = [new_end + n_start]

    return df


def reverse_df(df, n_start, n_stop, gene_index):
    """
    Given a neighborhood dataframe we want to reverse, performs the following:
    (i) Reverses orientation of each gene.
    (ii) Reverses order genes appear in (e.g. downstream genes [A, B, C] -> [C, B, A]).
    (iii) Reverses upstream/downstream genes (e.g. [upstream, AMR, downstream] -> [downstream, AMR, upstream]).
    (iv) Modify gene start and end positions to reflect neighborhood reversal.
    """
    # (iv) Reverse gene start/end positions to reflect reversal
    neighborhood_df = reverse_start_end(df, n_start, n_stop)

    # (i) Swap gene orientations
    swapped_df = swap_neighborhood_orientation(neighborhood_df)
    swapped_df.reset_index(drop=True, inplace=True)
    sorted_df = swapped_df.sort_values(by="Gene_Start")
    sorted_df.reset_index(drop=True, inplace=True)

    return swapped_df.copy(deep=True)


def make_neighborhood_JSON_data(gene_neighborhoods_dict, gene, neighborhood_size):
    """
    Creates dictionary of data required to write neighborhoods JSON file.
    This function should be run on the complete set of neighborhoods (i.e. for all AMR genes).
    """
    # Keeps track of data for clusters, links, and groups
    neighborhood_json_data = {}

    # CLUSTER DATA
    cluster_data = {}
    unique_genes = []

    # Genome gene contig dict
    genome_contigs = {}

    for genome_id, neighborhood_df in gene_neighborhoods_dict.items():
        neighborhood_data = {}
        contigs_dict = {}

        # cluster_uids: Assign a random alphanumeric string of 36 characters for uid (as per clustermap.js examples),
        neighborhood_data["uid"] = generate_alphanumeric_string(36)
        neighborhood_data["name"] = genome_id

        loci_data = {}
        loci_data["uid"] = generate_alphanumeric_string(36)
        loci_data["name"] = genome_id

        # Retain neighborhood start/end coords
        df = neighborhood_df.copy(deep=True)
        loci_data["start"] = df["Gene_Start"].min()
        loci_data["end"] = df["Gene_End"].max()

        # cluster_loci: Key: uid (gene index in neighborhood), Values: (arr) gene name, start coord, end coord, strand
        genes_dict = {}

        indices = neighborhood_df.index.values.tolist()

        # Reverse the neighborhood representation if the focal gene is not oriented correctly
        gene_index = int((len(indices) - 1) / 2)

        if neighborhood_df["Gene_Strand"][gene_index] == -1:
            neighborhood_df = reverse_df(
                neighborhood_df, loci_data["start"], loci_data["end"], gene_index
            )

        focal_gene_index = (neighborhood_size * 2 + 1) // 2
        for i in neighborhood_df.index:
            gene_data = {}
            gene_data["uid"] = neighborhood_df["Locus_Tag"][i].replace("'", "")
            if i == focal_gene_index:
                if gene not in unique_genes:
                    unique_genes.append(gene)
                gene_data["name"] = gene.replace('"', "")
            else:
                gene_data["name"] = neighborhood_df["Gene_Name"][i].replace('"', "")
            if gene_data["name"] not in unique_genes:
                unique_genes.append(gene_data["name"])
            gene_data["start"] = neighborhood_df["Gene_Start"][i]
            gene_data["end"] = neighborhood_df["Gene_End"][i]
            gene_data["strand"] = neighborhood_df["Gene_Strand"][i]

            # Keep track of contig data
            contigs_dict[neighborhood_df["Locus_Tag"][i]] = neighborhood_df[
                "Locus_Tag"
            ][i]

            genes_dict[i] = gene_data

        loci_data["genes"] = genes_dict
        neighborhood_data["loci"] = loci_data
        cluster_data[genome_id] = neighborhood_data
        genome_contigs[genome_id] = contigs_dict

    neighborhood_json_data["clusters"] = cluster_data

    # LINK DATA
    links_data = {}

    # GROUPS DATA
    groups_data = {}
    ind = 0
    for i in range(len(unique_genes)):
        gene = unique_genes[i]
        genomes = list(cluster_data.keys())

        # Ignore unidentified genes
        if not gene.startswith("UID"):
            group_data = {}
            gene_group_list = []

            for genome_1, genome_2 in itertools.combinations(genomes, 2):
                link_data = {}

                genome_1_presence = False
                genome_2_presence = False
                genome_1_key = 0
                genome_2_key = 0

                for id in cluster_data[genome_1]["loci"]["genes"].keys():
                    if gene in cluster_data[genome_1]["loci"]["genes"][id]["name"]:
                        genome_1_presence = True
                        genome_1_key = id
                        break

                for id in cluster_data[genome_2]["loci"]["genes"].keys():
                    if gene in cluster_data[genome_2]["loci"]["genes"][id]["name"]:
                        genome_2_presence = True
                        genome_2_key = id
                        break

                if genome_1_presence and genome_2_presence:
                    # Add link
                    target = {}
                    target["uid"] = cluster_data[genome_1]["loci"]["genes"][
                        genome_1_key
                    ]["uid"]
                    target["name"] = cluster_data[genome_1]["loci"]["genes"][
                        genome_1_key
                    ]["name"].replace("'", "")

                    query = {}
                    query["uid"] = cluster_data[genome_2]["loci"]["genes"][
                        genome_2_key
                    ]["uid"]
                    query["name"] = cluster_data[genome_2]["loci"]["genes"][
                        genome_2_key
                    ]["name"].replace("'", "")

                    contig_1 = target["uid"]
                    contig_2 = query["uid"]
                    link_data["uid"] = (
                        str(genome_1)
                        + "_"
                        + contig_1
                        + "-"
                        + str(genome_2)
                        + "_"
                        + contig_2
                    )
                    link_data["target"] = target
                    link_data["query"] = query
                    link_data["identity"] = 0.70

                    links_data[ind] = link_data
                    ind += 1

                if (
                    genome_1_presence
                    and not cluster_data[genome_1]["loci"]["genes"][genome_1_key]["uid"]
                    in gene_group_list
                ):
                    gene_group_list.append(
                        cluster_data[genome_1]["loci"]["genes"][genome_1_key]["uid"]
                    )
                if (
                    genome_2_presence
                    and not cluster_data[genome_2]["loci"]["genes"][genome_2_key]["uid"]
                    in gene_group_list
                ):
                    gene_group_list.append(
                        cluster_data[genome_2]["loci"]["genes"][genome_2_key]["uid"]
                    )

            group_data["uid"] = "group" + str(i)
            group_data["label"] = gene
            group_data["genes"] = gene_group_list
            groups_data[i] = group_data

    neighborhood_json_data["links"] = links_data
    neighborhood_json_data["groups"] = groups_data

    return neighborhood_json_data


def load_surrogates_file(output_path, gene):
    """
    Loads surrogates file into dictionary where keys are surrogate genome IDs, values are lists of genomes they are
    respectively standing in for.
    """
    surrogates_dict = {}
    with open(
        output_path + "/JSON/surrogates/" + gene + "_surrogates.txt", "r"
    ) as infile:
        data = infile.readlines()

    for line in data:
        # Case 1: Surrogate genome with list of represented genomes
        if ":" in line:
            surrogates_data = line.split(": ")
            surrogates = surrogates_data[1]
            surrogates_list = (
                surrogates.strip().replace("[", "").replace("]", "").split(", ")
            )
            surrogates_dict[surrogates_data[0]] = surrogates_list

        # Case 2: Singleton representative genome
        else:
            surrogates_dict[line.strip("\n")] = []

    return surrogates_dict


def write_neighborhood_JSON(
    neighborhood_json_data, gene, output_path, surrogates=False
):
    """
    Creates JSON file containing neighborhood data for an AMR gene for the genomes under analysis of the following form
    (where the first N neighbor genes represent the N neighbors downstream of the target gene, and the last N
    neighbor genes represent the N neighbors upstream of the target gene).
    """
    # Write JSON format
    out_path = output_path + "/JSON"
    check_output_path(out_path)

    try:
        cluster_data = neighborhood_json_data["clusters"]
        links_data = neighborhood_json_data["links"]
        groups_data = neighborhood_json_data["groups"]
    except KeyError:
        print("Neighborhood data for {g} was not found.".format(g=gene))
        return

    if surrogates:
        output_file_path = out_path + "/" + gene + "_surrogates.json"
    else:
        output_file_path = out_path + "/" + gene + "_temp.json"

    with open(output_file_path, "w") as outfile:
        outfile.write("{\n")

        # ----------------------------------------------- CLUSTERS -----------------------------------------------------
        outfile.write('\t"clusters": [\n')

        cluster_index = 0
        num_clusters = len(cluster_data.keys())
        for genome_id in cluster_data.keys():
            outfile.write("\t\t{\n")

            # Neighborhood uid and name
            outfile.write(
                '\t\t\t"uid": ' + '"' + cluster_data[genome_id]["uid"] + '",\n'
            )
            outfile.write('\t\t\t"name": ' + '"' + genome_id + '",\n')

            # Neighborhood details: uid, name, start, stop, end, and gene data
            outfile.write('\t\t\t"loci": [\n')
            outfile.write("\t\t\t\t{\n")
            outfile.write(
                '\t\t\t\t\t"uid": '
                + '"'
                + cluster_data[genome_id]["loci"]["uid"]
                + '",\n'
            )
            outfile.write(
                '\t\t\t\t\t"name": '
                + '"'
                + cluster_data[genome_id]["loci"]["name"]
                + '",\n'
            )
            outfile.write(
                '\t\t\t\t\t"start": '
                + str(cluster_data[genome_id]["loci"]["start"])
                + ",\n"
            )
            outfile.write(
                '\t\t\t\t\t"end": '
                + str(cluster_data[genome_id]["loci"]["end"])
                + ",\n"
            )

            # Neighborhood genes data
            outfile.write('\t\t\t\t\t"genes": [\n')

            index = 0
            neighborhood_length = len(cluster_data[genome_id]["loci"]["genes"])
            for gene_key, gene_data in cluster_data[genome_id]["loci"]["genes"].items():
                # Opening bracket
                outfile.write("\t\t\t\t\t\t{\n")

                # Gene data
                outfile.write(
                    '\t\t\t\t\t\t\t"uid": '
                    + '"'
                    + cluster_data[genome_id]["loci"]["genes"][gene_key]["uid"]
                    + '",\n'
                )
                outfile.write(
                    '\t\t\t\t\t\t\t"name": '
                    + '"'
                    + cluster_data[genome_id]["loci"]["genes"][gene_key][
                        "name"
                    ].replace("'", "")
                    + '",\n'
                )
                outfile.write(
                    '\t\t\t\t\t\t\t"start": '
                    + str(cluster_data[genome_id]["loci"]["genes"][gene_key]["start"])
                    + ",\n"
                )
                outfile.write(
                    '\t\t\t\t\t\t\t"end": '
                    + str(cluster_data[genome_id]["loci"]["genes"][gene_key]["end"])
                    + ",\n"
                )
                outfile.write(
                    '\t\t\t\t\t\t\t"strand": '
                    + str(cluster_data[genome_id]["loci"]["genes"][gene_key]["strand"])
                    + "\n"
                )

                # Closing bracket
                if index != len(cluster_data[genome_id]["loci"]["genes"]) - 1:
                    outfile.write("\t\t\t\t\t\t},\n")
                    index += 1
                else:
                    outfile.write("\t\t\t\t\t\t}\n")

            # Closing brackets for clusters
            outfile.write("\t\t\t\t\t]\n")
            outfile.write("\t\t\t\t}\n")
            outfile.write("\t\t\t]\n")
            if cluster_index != num_clusters - 1:
                outfile.write("\t\t},\n")
                cluster_index += 1
            else:
                outfile.write("\t\t}\n")

        # Clusters data final closing bracket
        outfile.write("\t],\n")

        # ----------------------------------------------- LINKS --------------------------------------------------------
        outfile.write('\t"links": [\n')

        link_index = 0
        num_links = len(links_data.keys())
        last_link_flag = False
        for link in links_data.keys():
            try:
                uid = links_data[link]["target"]["uid"]
                outfile.write("\t\t{\n")
                outfile.write('\t\t\t"uid": "' + links_data[link]["uid"] + '",\n')
                outfile.write('\t\t\t"target": {\n')
                outfile.write(
                    '\t\t\t\t"uid": "' + links_data[link]["target"]["uid"] + '",\n'
                )
                outfile.write(
                    '\t\t\t\t"name": "' + links_data[link]["target"]["name"] + '-1"\n'
                )
                outfile.write("\t\t\t},\n")
                outfile.write('\t\t\t"query": {\n')
                outfile.write(
                    '\t\t\t\t"uid": "' + links_data[link]["query"]["uid"] + '",\n'
                )
                outfile.write(
                    '\t\t\t\t"name": "' + links_data[link]["query"]["name"] + '-2"\n'
                )
                outfile.write("\t\t\t},\n")
                outfile.write(
                    '\t\t\t"identity": "' + str(links_data[link]["identity"]) + '"\n'
                )

            except KeyError:
                pass

            if link_index != num_links - 1:
                outfile.write("\t\t},\n")

            else:
                outfile.write("\t\t}\n")
                last_link_flag = True

            link_index += 1

            if last_link_flag:
                break

        # Links data final closing bracket
        outfile.write("\t],\n")

        # ------------------------------------------------ GROUPS ------------------------------------------------------
        outfile.write('\t"groups": [\n')

        group_index = 0
        num_groups = len(groups_data.keys())
        for group in groups_data.keys():
            outfile.write("\t\t{\n")

            outfile.write('\t\t\t"uid": "' + groups_data[group]["uid"] + '",\n')
            outfile.write('\t\t\t"label": "' + groups_data[group]["label"] + '",\n')
            outfile.write(
                '\t\t\t"genes": ' + json.dumps(groups_data[group]["genes"]) + "\n"
            )

            if group_index != num_groups - 1:
                outfile.write("\t\t},\n")
                group_index += 1
            else:
                outfile.write("\t\t}\n")

        outfile.write("\t]\n")

        # Final closing bracket
        outfile.write("}\n")


def write_clustermap_JSON_HTML(gene, sample_data_path, out_path, rep_type="standard"):
    """
    Generates accompanying HTML file for clustermap compatible JSON representation of neighborhood.
    Creates standalone HTML file for each respective type of neighborhood representation (e.g. standard,
    surrogates, or with representative UPGMA cluster) in case user wants to load individual visualizations.
    """
    json_filename = ""
    html_filename = ""

    if rep_type == "upgma":
        json_filename = f"{gene}_upgma.json"
        html_filename = f"{gene}_upgma.html"
    elif rep_type == "surrogates":
        json_filename = f"{gene}_surrogates.json"
        html_filename = f"{gene}_surrogates.html"
    else:  # rep_type == 'standard'
        json_filename = f"{gene}_temp.json"
        html_filename = f"{gene}.html"

    json_path = f"{out_path}/JSON/{json_filename}"
    file_path = f"{out_path}/JSON/{html_filename}"
    second_line = f'\t\td3.json("{json_filename}")\n'

    # Make HTML index file with appropriate JSON
    with open(file_path, "w") as html_outfile, open(sample_data_path) as template:
        for line in template:
            if "height: 100vh;" in line:
                # Determine number of genomes present
                with open(json_path, "r") as infile:
                    if len(infile.readlines()) != 0:
                        infile.seek(0)
                        json_data = json.load(infile)
                num_clusters = len(json_data["clusters"])

                # Calculate optimal canvas height as proportional to number of genomes
                height = int(64.67 * num_clusters + 100)
                height_px = str(height) + "px"
                if height < 900:
                    height_px = "100vh"
                html_outfile.write("\t\t\t\theight: {};\n".format(height_px))
            else:
                html_outfile.write(line)

        html_outfile.write("\n")
        html_outfile.write(second_line)
        html_outfile.write("\t\t\t.then(data => {\n")
        html_outfile.write('\t\t\t\tdiv.selectAll("div")\n')
        html_outfile.write("\t\t\t\t\t.data([data])\n")
        html_outfile.write('\t\t\t\t\t.join("div")\n')
        html_outfile.write("\t\t\t\t\t.call(chart)\n\n")
        html_outfile.write('\t\t\t\tlet svg = div.select("svg")\n')
        html_outfile.write('\t\t\t\td3.select("#btn-save-svg")\n')
        html_outfile.write('\t\t\t\t\t.on("click", () => {\n')
        html_outfile.write("\t\t\t\t\t\tconst blob = serialise(svg)\n")
        html_outfile.write('\t\t\t\t\t\tdownload(blob, "clinker.svg")\n')
        html_outfile.write("\t\t\t\t\t})\n")
        html_outfile.write("\t\t\t})\n")
        html_outfile.write("\t</script>\n")
        html_outfile.write("</html>")


def make_gene_HTML(genes_list, sample_data_path, out_path):
    """
    For each AMR gene for which a JSON file was created, generates an accompanying HTML file for rendering its gene
    order visualization using clustermap with. This is done for each gene individually.
    """
    for gene in genes_list:
        # Make HTML index file with appropriate JSON
        write_clustermap_JSON_HTML(gene, sample_data_path, out_path)


def get_cluster_data_genes_uid_list(json_cluster_data, genome_ids):
    """
    Given cluster data for a JSON clustermap representation of a gene's neighborhoods and a list of genomes to extract
    from, returns a list containing every gene uid present in the neighborhoods.
    """
    gene_uids_list = []
    for genome_id in genome_ids:
        for gene_key, gene_data in json_cluster_data[genome_id]["loci"][
            "genes"
        ].items():
            gene_uids_list.append(
                json_cluster_data[genome_id]["loci"]["genes"][gene_key]["uid"]
            )
    return gene_uids_list


def remove_defunct_clustermap_data(json_data):
    """
    For UPGMA representative clustermap JSON representations created based on original JSON file, after removing all
    cluster data for genomes no longer included, this function allows us to delete residual links and gene groups that
    no longer apply.
    """
    genome_ids = []
    for cluster in json_data["clusters"]:
        genome_ids.append(cluster["name"])
    gene_links = []
    for link in json_data["links"]:
        # Get both genome names
        genome_contig_details = link["uid"].split("-")
        genome_1 = genome_contig_details[0].split("_")[0]
        genome_2 = genome_contig_details[1].split("_")[0]

        if genome_1 in genome_ids and genome_2 in genome_ids:
            gene_links.append(link)

    gene_uids_list = get_cluster_data_genes_uid_list(json_data["clusters"], genome_ids)
    gene_groups = []
    for group in json_data["groups"]:
        if group["uid"] in gene_uids_list:
            gene_groups.append(group)

    json_data["links"] = gene_links
    json_data["groups"] = gene_groups

    return json_data


def rename_clustermap_genes_contigs(output_path, fasta_path, gene, neighborhood_size):
    """
    Given clustermap-style JSON representing neighborhood data for a given gene across multiple genomes, modifies JSON
    to rename all gene uids according to their contig. This is required for downstream filtering to obtain unique
    neighborhoods based on percent identity of matched hits, wherein identical neighborhoods are represented by a
    surrogate neighborhood.
    """
    # Load gene JSON surrogate representation
    json_data, gene_path = load_JSON_data(output_path, gene, json_file="surrogates")

    # Get contig dictionary for gene
    fasta_contigs_dict = make_fasta_contig_dict(fasta_path, gene)

    # Load neighborhood indices
    neighborhood_indices_dict = json.load(
        open(output_path + "/neighborhood_indices.txt")
    )

    for cluster in json_data["clusters"]:
        # Loop over each locus in the cluster
        for locus in cluster["loci"]:
            # Get neighborhood indices and adjust them
            genomes_indices = neighborhood_indices_dict[gene]
            neighborhood_indices = genomes_indices[cluster["name"]]

            # Case 1: no upstream genes
            if neighborhood_indices[0] == 0:
                focal_gene_index = 0
            # Case 2: no downstream genes
            elif neighborhood_indices[1] == 0:
                focal_gene_index = 0
            # Case 3: full neighborhood
            elif (
                neighborhood_indices[0] == -neighborhood_size
                and neighborhood_indices[1] == neighborhood_size
            ):
                focal_gene_index = (neighborhood_indices[1] * 2) // 2
            # Case 4: at least one upstream, downstream gene
            else:
                neighborhood_indices[0] += abs(neighborhood_indices[0])
                neighborhood_indices[1] += abs(neighborhood_indices[0])
                focal_gene_index = neighborhood_indices[1] // 2

            # Get the orientation of the focal gene
            orientation = locus["genes"][focal_gene_index]["strand"]

            # Loop over each gene in the locus
            for gene_entry in locus["genes"]:
                # Get the gene index from its uid and convert to an integer
                gene_uid_tokens = gene_entry["uid"].split("-")
                gene_index = int(gene_uid_tokens[1])

                # Get the contig ID for this gene from the fasta_contigs_dict
                contig_id = fasta_contigs_dict[cluster["name"]][gene_index]

                # Replace the gene uid with the contig ID
                gene_entry["uid"] = contig_id

    # Save new data
    with open(gene_path, "w") as outfile:
        json.dump(json_data, outfile)


def load_JSON_data(output_path, gene, json_file="base"):
    """
    Helper function for loading JSON neighborhood data.
    """
    json_data = ""
    if json_file == "surrogates":
        gene_path = output_path + "/JSON/" + gene + "_surrogates.json"
    elif json_file == "temp":
        gene_path = output_path + "/JSON/" + gene + "_temp.json"
    else:
        gene_path = output_path + "/JSON/" + gene + ".json"

    with open(gene_path, "r") as infile:
        if len(infile.readlines()) != 0:
            infile.seek(0)
            json_data = json.load(infile)

    return json_data, gene_path


def update_JSON_links_PI(BLAST_df_dict, output_path):
    """
    Updates JSON representations of AMR gene neighborhoods created using extraction module so that gene cluster links
    reflect percent identities found in blast results.
    """
    for gene, blast_files_dict in BLAST_df_dict.items():
        # Load AMR gene JSON link data
        json_data, original_gene_path = load_JSON_data(
            output_path, gene, json_file="temp"
        )
        original_gene_path = Path(original_gene_path)

        final_gene_path = output_path + "/JSON/" + gene + ".json"

        # Update each link according to the respective blast results
        for i in range(len(json_data["links"])):
            # Contig identifiers
            locus_data = json_data["links"][i]["uid"].split("-")
            locus_1 = locus_data[0].split("_")
            locus_2 = locus_data[1].split("_")
            locus_tag_1 = locus_1[1] + "_" + locus_1[2]
            locus_tag_2 = locus_2[1] + "_" + locus_2[2]

            # Genome names
            genome_1 = locus_1[0]
            genome_2 = locus_2[0]

            # Check respective BLAST file dataframe and update percent identity
            if genome_1 + "_" + genome_2 + ".blast.txt" in blast_files_dict.keys():
                df = blast_files_dict[genome_1 + "_" + genome_2 + ".blast.txt"]
                row = df.filter(
                    ((df["query_id"] == locus_tag_1) & (df["sub_id"] == locus_tag_2))
                    | (df["query_id"] == locus_tag_2) & (df["sub_id"] == locus_tag_1)
                )
            else:
                df = blast_files_dict[genome_2 + "_" + genome_1 + ".blast.txt"]
                row = df.filter(
                    ((df["query_id"] == locus_tag_1) & (df["sub_id"] == locus_tag_2))
                    | (df["query_id"] == locus_tag_2) & (df["sub_id"] == locus_tag_1)
                )

            try:
                PI = row.get_column("PI").item()
            except ValueError:
                PI = 0.70

            json_data["links"][i]["identity"] = PI

        # Overwrite JSON file with updated data
        with open(final_gene_path, "w") as outfile:
            json.dump(json_data, outfile)

        # Remove 'temp' JSON file
        original_gene_path.unlink(missing_ok=True)


def map_genome_id_to_dendrogram_leaves(upgma_clusters, genome_to_num_mapping):
    """
    Given a UPGMA dendrogram and a dictionary mapping genomes to integers corresponding to the indices of the distance
    matrix used to generate the dendrogram, returns a list of genome names in order of the dendrogram leaves from left
    to right.
    """
    upgma_leaves = upgma_clusters["ivl"]
    genome_order = []
    for leaf in upgma_leaves:
        genome = genome_to_num_mapping[leaf]
        genome_order.append(genome)

    return genome_order


def order_cluster_data_by_dendrogram(genome_order_dict, json_cluster_data):
    """
    Given JSON cluster data in the clustermap format, reorders cluster data according to dendrogram leaves from
    left to right.
    """
    clusters = []
    for genome in genome_order_dict:
        for cluster in json_cluster_data:
            if cluster["name"] == genome:
                clusters.append(cluster)
    return clusters


def order_JSON_clusters_UPGMA(
    output_path, gene, upgma_clusters, genome_to_num_mapping, surrogates=False
):
    """
    Reorders how genomes are encoded in their respective JSON files for an AMR gene according to how they
    were clustered by UPGMA (from left to right).
    """
    # Load AMR gene JSON cluster data
    json_data, gene_path = load_JSON_data(output_path, gene, surrogates)

    # Reorder cluster data genomes according to UPGMA leaves ordering
    genome_order = map_genome_id_to_dendrogram_leaves(
        upgma_clusters, genome_to_num_mapping
    )
    clusters = order_cluster_data_by_dendrogram(genome_order, json_data["clusters"])

    # Update JSON data
    json_data["clusters"] = clusters

    # Update file
    with open(gene_path, "w") as outfile:
        json.dump(json_data, outfile)


def update_cluster_data(json_data, representative_cluster_genomes):
    """
    Given some JSON data and a list of representative genomes, removes data for all other genomes.
    """
    clusters = []
    for genome in representative_cluster_genomes:
        for cluster in json_data["clusters"]:
            if cluster["name"] == genome:
                clusters.append(cluster)

    # Update JSON data
    json_data["clusters"] = clusters

    return json_data


def clean_json_data(json_data):
    """
    For UPGMA representative clustermap JSON representations created based on original JSON file, after removing all
    cluster data for genomes no longer included, this function allows us to delete residual links and gene groups that
    no longer apply.
    """
    genome_ids = []
    for cluster in json_data["clusters"]:
        genome_ids.append(cluster["name"])
    gene_links = []
    for link in json_data["links"]:
        # Get both genome names
        genome_contig_details = link["uid"].split("-")
        genome_1 = genome_contig_details[0].split("_")[0]
        genome_2 = genome_contig_details[1].split("_")[0]

        if genome_1 in genome_ids and genome_2 in genome_ids:
            gene_links.append(link)

    json_data["links"] = gene_links

    return json_data


def make_representative_UPGMA_cluster_JSON(
    output_path, gene, upgma_clusters, genome_to_num_mapping
):
    """
    Creates a JSON file with one representative genome from each UPGMA cluster.
    """
    # Load AMR gene JSON cluster data
    json_data = ""
    with open(output_path + "/JSON/" + gene + ".json", "r") as infile:
        if len(infile.readlines()) != 0:
            infile.seek(0)
            json_data = json.load(infile)

    # Determine genomes to keep in representation
    representative_cluster_genomes = []

    genome_order = map_genome_id_to_dendrogram_leaves(
        upgma_clusters, genome_to_num_mapping
    )
    cluster_df = pd.DataFrame(
        {"genome": genome_order, "cluster": upgma_clusters["leaves_color_list"]}
    )

    unique_clusters = set(upgma_clusters["leaves_color_list"])
    for cluster in unique_clusters:
        cluster_rows = cluster_df.loc[cluster_df["cluster"] == cluster]
        representative_cluster_genomes.append(cluster_rows.iloc[0].genome)
    del cluster_df

    # Update cluster data to only include those genomes
    json_data = update_cluster_data(json_data, representative_cluster_genomes)

    # Get rid of defunct links and gene groups
    final_json_data = clean_json_data(json_data)

    # Update file
    with open(output_path + "/JSON/" + gene + "_upgma.json", "w+") as outfile:
        json.dump(final_json_data, outfile)

    # Make respective HTML file for Coeus
    write_clustermap_JSON_HTML(
        gene, output_path + "/index.html", output_path, rep_type="upgma"
    )
