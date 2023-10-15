#!/usr/bin/env python
"""
Methods for generating JSON representations of extracted neighborhoods that can be used with clustermap.js.
"""
import os
import re
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


def get_focal_index(neighborhood_indices, num_neighbors):
    """
    Given a genome's neighborhood indices list as retrieved from the neighborhood_indices.txt output from extraction.py,
    determines the index of the focal gene.
    """
    # Case 1: no upstream genes
    if neighborhood_indices[0] == 0:
        focal_gene_index = 0

    # Case 2: no downstream genes
    elif neighborhood_indices[1] == 0:
        focal_gene_index = 0

    # Case 3: full neighborhood
    elif (
        neighborhood_indices[0] == -num_neighbors
        and neighborhood_indices[1] == num_neighbors
    ):
        focal_gene_index = (neighborhood_indices[1] * 2) // 2

    # Case 4: at least one upstream, downstream gene
    else:
        focal_gene_index = abs(neighborhood_indices[0])

    return focal_gene_index


def make_neighborhood_JSON_data(
    gene_neighborhoods_dict, neighborhood_indices, gene, neighborhood_size
):
    """
    Creates dictionary of data required to write neighborhoods JSON file.
    This function should be run on the complete set of neighborhoods (i.e. for all AMR genes).
    """
    # Keeps track of data for clusters, links, and groups
    neighborhood_json_data = {}

    # CLUSTER DATA
    cluster_data = {}
    reversed_genomes = {}
    unique_genes = []

    # Genome gene contig dict
    genome_contigs = {}

    for genome_id, neighborhood_df in gene_neighborhoods_dict.items():
        neighborhood_data = {}
        contigs_dict = {}

        neighborhood_data["uid"] = generate_alphanumeric_string(36)
        neighborhood_data["name"] = genome_id

        loci_data = {}
        loci_data["uid"] = generate_alphanumeric_string(36)
        loci_data["name"] = genome_id

        # Retain neighborhood start/end coords
        df = neighborhood_df.copy(deep=True)
        loci_data["start"] = df["Gene_Start"].min()
        loci_data["end"] = df["Gene_End"].max()

        genes_dict = {}

        indices = neighborhood_df.index.values.tolist()

        # Reverse the neighborhood representation if the focal gene is not oriented correctly
        indices = neighborhood_indices[gene][genome_id]
        gene_index = get_focal_index(indices, neighborhood_size)

        if neighborhood_df["Gene_Strand"][gene_index] == -1:
            neighborhood_df = reverse_df(
                neighborhood_df, loci_data["start"], loci_data["end"], gene_index
            )
            indices = [-1 * indices[1], indices[0]]
            reversed_genomes[genome_id] = True
        else:
            reversed_genomes[genome_id] = False

        focal_gene_index = get_focal_index(indices, neighborhood_size)
        neighborhood_df.reset_index(drop=True, inplace=True)

        for i in neighborhood_df.index:
            gene_data = {}
            gene_data["uid"] = neighborhood_df["Locus_Tag"][i].replace("'", "")
            gene_data["name"] = neighborhood_df["Gene_Name"][i].replace('"', "")

            # Add gene name to list for group names: if focal gene, ignore any extra annotation details
            if i == focal_gene_index:
                if gene not in unique_genes:
                    unique_genes.append(gene)
            else:
                if (
                    gene not in gene_data["name"]
                    and gene_data["name"] not in unique_genes
                ):
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

    ind = 0
    for i in range(len(unique_genes)):
        unique_gene = unique_genes[i]
        genomes = list(cluster_data.keys())

        # Ignore unidentified genes
        if not unique_gene.startswith("UID"):
            for genome_1, genome_2 in itertools.combinations(genomes, 2):
                link_data = {}

                genome_1_presence = False
                genome_2_presence = False
                genome_1_key = 0
                genome_2_key = 0

                for id in cluster_data[genome_1]["loci"]["genes"].keys():
                    if (
                        unique_gene
                        in cluster_data[genome_1]["loci"]["genes"][id]["name"]
                    ):
                        genome_1_presence = True
                        genome_1_key = id
                        break

                for id in cluster_data[genome_2]["loci"]["genes"].keys():
                    if (
                        unique_gene
                        in cluster_data[genome_2]["loci"]["genes"][id]["name"]
                    ):
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
                    link_data["uid"] = f"{genome_1}_{contig_1}|{genome_2}_{contig_2}"
                    link_data["target"] = target
                    link_data["query"] = query
                    if gene in target["name"] and gene in query["name"]:
                        link_data["identity"] = 0
                    else:
                        link_data["identity"] = 0.70

                    links_data[ind] = link_data
                    ind += 1

    # GROUPS DATA
    groups_data = {}
    group_id = 1

    for i in range(len(unique_genes)):
        unique_gene = unique_genes[i]
        genomes = list(cluster_data.keys())

        # Ignore unidentified genes
        if not unique_gene.startswith("UID"):
            group_data = {}
            gene_uids = set()
            for genome in genomes:
                for id in cluster_data[genome]["loci"]["genes"].keys():
                    if unique_gene in cluster_data[genome]["loci"]["genes"][id]["name"]:
                        gene_uids.add(cluster_data[genome]["loci"]["genes"][id]["uid"])

            group_data["uid"] = f"group{group_id}"
            group_data["label"] = unique_gene
            group_data["genes"] = list(gene_uids)

            groups_data[unique_gene] = group_data
            group_id += 1

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


def get_JSON_specific_configs(json_filename, json_path):
    """
    Determines values for json filepath and clustermap chart height to write to JSON HTML file.
    """
    json_file_path = f'\t\td3.json("{json_filename}")\n'

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

    return json_file_path, height_px


def write_clustermap_JSON_HTML(gene, out_path, rep_type="standard"):
    """
    Generates accompanying HTML file for clustermap compatible JSON representation of neighborhood.
    Creates standalone HTML file for each respective type of neighborhood representation (e.g. standard,
    surrogates, or with representative UPGMA cluster) in case user wants to load individual visualizations.
    """

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

    json_file_path, height_px = get_JSON_specific_configs(json_filename, json_path)

    # Write HTML contents to file to represent clustermap chart for the gene
    html_content = """\
    <!DOCTYPE html>
    <html>
        <head>
            <meta charset="utf-8">
            <title>cmap</title>
            <script src="../dist/d3.min.js"></script>
            <style>
                body {{ margin: 0; padding: 0; }}
                div {{
                    width: 100vw;
                    height: {height_value};
                    margin: 0;
                    padding: 0;
                }}
            </style>
        </head>
        <body>
            <main>
                <button id="btn-save-svg">Save</button>
                <div id="plot"></div>
            </main>
        </body>
        <script type="module">
            import clusterMap from "../src/clusterMap.js"
            function serialise(svg) {{
                /* Saves the figure to SVG in its current state.
                 * Clones the provided SVG and sets the width/height of the clone to the
                 * bounding box of the original SVG. Thus, downloaded figures will be sized
                 * correctly.
                 * This function returns a new Blob, which can then be downloaded.
                */
                let node = svg.node();
                const xmlns = "http://www.w3.org/2000/xmlns/";
                const xlinkns = "http://www.w3.org/1999/xlink";
                const svgns = "http://www.w3.org/2000/node";
                const bbox = svg.select("g").node().getBBox()
                node = node.cloneNode(true);
                node.setAttribute("width", bbox.width);
                node.setAttribute("height", bbox.height);
                node.setAttributeNS(xmlns, "xmlns", svgns);
                node.setAttributeNS(xmlns, "xmlns:xlink", xlinkns);
                // Adjust x/y of <g> to account for axis/title position.
                // Replaces the transform attribute, so drag/zoom is ignored.
                d3.select(node)
                    .select("g")
                    .attr("transform", `translate({{Math.abs(bbox.x)}}, {{Math.abs(bbox.y)}})`)
                const serializer = new window.XMLSerializer;
                const string = serializer.serializeToString(node);
                return new Blob([string], {{type: "image/node+xml"}});
            }}
            function download(blob, filename) {{
                /* Downloads a given blob to filename.
                 * This function appends a new anchor to the document, which points to the
                 * supplied blob. The anchor.click() method is called to trigger the download,
                 * then the anchor is removed.
                */
                const link = document.createElement("a");
                link.href = URL.createObjectURL(blob);
                link.download = filename;
                document.body.appendChild(link);
                link.click();
                document.body.removeChild(link);
            }}
            const div = d3.select("#plot")
                .attr("width", "2400vw")
                .attr("height", "{height_value}")
            const chart = clusterMap()
                .config({{
                    cluster: {{
                        alignLabels: true
                    }},
                    gene: {{
                        label: {{
                            show: false,
                        }}
                    }},
                    link: {{
                        threshold: 0.3,
                        bestOnly: true,
                    }}
                }})
            d3.json("{path_to_json}")
                    .then(data => {{
                        div.selectAll("div")
                            .data([data])
                            .join("div")
                            .call(chart)
                        let svg = div.select("svg")
                        d3.select("#btn-save-svg")
                            .on("click", () => {{
                                const blob = serialise(svg)
                                download(blob, "clinker.svg")
                            }})
                    }})
        </script>
    </html>
    """

    with open(file_path, "w") as html_outfile:
        html_outfile.write(
            html_content.format(height_value=height_px, path_to_json=json_file_path)
        )


def make_gene_HTML(genes_list, out_path):
    """
    For each AMR gene for which a JSON file was created, generates an accompanying HTML file for rendering its gene
    order visualization using clustermap with. This is done for each gene individually.
    """
    for gene in genes_list:
        # Make HTML index file with appropriate JSON
        write_clustermap_JSON_HTML(gene, out_path)


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


def add_UID_genes_links_groups(BLAST_df_dict, output_path):
    """
    Adds links between all unidentified genes. Initially sets PI to 0, as we update the PI (if found) at a later
    step and otherwise delete the link.
    """
    json_dict = dict()
    for gene, blast_files_dict in BLAST_df_dict.items():
        # Load AMR gene JSON link data
        json_data, gene_path = load_JSON_data(output_path, gene, json_file="temp")

        # Fetch link data and determine number of existing links for counter
        clusters = json_data["clusters"]
        links = json_data["links"]

        # Create a dictionary to store genes by name for each cluster
        genes_by_name = {}

        # Iterate over clusters and populate the dictionary
        for cluster in clusters:
            genes_by_name[cluster["name"]] = [
                gene
                for gene in cluster["loci"][0]["genes"]
                if gene["name"].startswith("UID")
            ]

        # Iterate over every combination of clusters
        for cluster_1, cluster_2 in itertools.combinations(clusters, 2):
            genes_1 = genes_by_name[cluster_1["name"]]
            genes_2 = genes_by_name[cluster_2["name"]]

            # Iterate over unidentified genes and add temporary link between all unique combinations
            for gene_1 in genes_1:
                for gene_2 in genes_2:
                    link_uid = f"{cluster_1['name']}_{gene_1['uid']}|{cluster_2['name']}_{gene_2['uid']}"

                    # Create a new link object
                    link = {
                        "uid": link_uid,
                        "target": {"uid": gene_1["uid"], "name": gene_1["name"]},
                        "query": {"uid": gene_2["uid"], "name": gene_2["name"]},
                        "identity": "0",
                    }

                    # Add the link to the links list
                    links.append(link)

        # Update links
        json_data["links"] = links
        json_dict[gene] = json_data

    return json_dict


def update_UID_groups(json_data):
    """
    Adds groups for unidentified genes so they can be toggled in clustermap legend.
    """
    # Create a list of unidentified genes across all clusters
    unidentified_genes_by_genome = dict()
    unidentified_genes = []
    for cluster in json_data["clusters"]:
        for gene in cluster["loci"][0]["genes"]:
            if gene["name"].startswith("UID"):
                unidentified_genes.append(gene)
                unidentified_genes_by_genome[gene["uid"]] = cluster["name"]

    # Create a shallow copy, candidate_genes, to avoid iterating over genes already added to a group
    candidate_genes = unidentified_genes[:]

    group_index = len(json_data["groups"]) + 1
    uid_group_index = 1
    for unidentified_gene in unidentified_genes:
        if unidentified_gene in candidate_genes:
            # Keep track of UID genes that are 70% or more identical
            group_genes = []
            unidentified_gene_uid = unidentified_gene["uid"]

            # Check links for percent identity >= 0.70
            for link in json_data["links"]:
                query_uid = link["query"]["uid"]
                target_uid = link["target"]["uid"]

                query_name = link["query"]["name"]
                target_name = link["target"]["name"]

                if (
                    unidentified_gene_uid == query_uid
                    and target_name.startswith("UID")
                    and unidentified_genes_by_genome[query_uid]
                    != unidentified_genes_by_genome[target_uid]
                ) or (
                    unidentified_gene_uid == target_uid
                    and query_name.startswith("UID")
                    and unidentified_genes_by_genome[target_uid]
                    != unidentified_genes_by_genome[query_uid]
                ):
                    # Add both genes to the new UID group if similar enough
                    if link["identity"] >= 0.70:
                        group_genes.append(target_uid)
                        group_genes.append(query_uid)

                        # Remove genes from candidate genes list
                        if query_name in candidate_genes:
                            candidate_genes.remove(query_name)
                        if target_name in candidate_genes:
                            candidate_genes.remove(target_name)

            # Create a new group if group_genes is not empty
            if group_genes:
                group = {
                    "uid": f"group{group_index}",
                    "label": f"UID Group {uid_group_index}",
                    "genes": list(set(group_genes)),
                }
                json_data["groups"].append(group)
                group_index += 1
                uid_group_index += 1

    return json_data


def make_assembly_dict(assembly_path):
    """
    Makes dictionary consisting of assembly file names with keys as the filename without the file extension, values
    as the full filename.
    """
    faa_dict = {}
    for filename in assembly_path:
        key = filename.split(".")[0]
        faa_dict[key] = filename

    return faa_dict


def update_JSON_links_PI(BLAST_df_dict, json_dict, output_path):
    """
    Updates JSON representations of AMR gene neighborhoods created using extraction module so that gene cluster links
    reflect percent identities found in blast results.
    """
    for gene, blast_files_dict in BLAST_df_dict.items():
        # Load AMR gene JSON link data
        old_json_data, original_gene_path = load_JSON_data(
            output_path, gene, json_file="temp"
        )
        json_data = json_dict[gene]
        original_gene_path = Path(original_gene_path)

        final_gene_path = output_path + "/JSON/" + gene + ".json"

        # Update each link according to the respective blast results
        updated_links = []
        for i in range(len(json_data["links"])):
            val_error_flag = False
            # Contig identifiers
            locus_data = json_data["links"][i]["uid"].split("|")
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
                val_error_flag = True

            # Delete links between annotated focal genes that are not homologous
            if not (json_data["links"][i]["identity"] == "0" and val_error_flag):
                json_data["links"][i]["identity"] = PI
                updated_links.append(json_data["links"][i])

        json_data["links"] = updated_links

        # Update UID group
        json_data = update_UID_groups(json_data)

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


def order_JSON_clusters_UPGMA(output_path, gene, upgma_clusters, genome_to_num_mapping):
    """
    Reorders how genomes are encoded in their respective JSON files for an AMR gene according to how they
    were clustered by UPGMA (from left to right).
    """
    # Load AMR gene JSON cluster data
    json_data, gene_path = load_JSON_data(output_path, gene)

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
        genome_contig_details = link["uid"].split("|")
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
    write_clustermap_JSON_HTML(gene, output_path, rep_type="upgma")
