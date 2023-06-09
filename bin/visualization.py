#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import random

from scipy.cluster import hierarchy
from skbio.stats.ordination import pcoa

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors

import plotly
import plotly.express as px
import plotly.graph_objects as go
from plotly.figure_factory import create_dendrogram

import networkx as nx
from matplotlib.pylab import savefig, cm, axis

def get_colors(num):
    """
    Color mapping for visualizations.
    """
    colors = ['#6e40aa', '#b83cb0', '#c33dad', '#ff4f7c', '#f6478d', '#ff6956', '#f59f30', '#c4d93e',
              '#83f557', '#38f17a', '#22e599', '#19d3b5', '#29a0dd', '#5069d9', '#5f56c9', '#bbbbbb']

    if num <= len(colors):
        return colors
    else:
        for i in range(num - len(colors)):
            rand_color = lambda: random.randint(0, 255)
            colors.append('#%02X%02X%02X' % (rand_color(), rand_color(), rand_color()))
        return colors


def plot_similarity_histogram(similarity_dict, save_path):
    """
    For use with clustering and scoring functionalities to plot a histogram of
    either average similarity scores across all extracted AMR gene neighborhoods.
    """
    similarity_scores = list(similarity_dict.values())

    # Number of items equivalent to number of AMR gene neighborhoods extracted
    num_points = len(similarity_dict)

    # Similarity and distance values will always be between 0 and 1: divide into 20 bins
    num_bins = int(round(max(list(similarity_dict.values()))))

    # Defining histogram
    fig, ax = plt.subplots(1, 1)
    N, bins, patches = ax.hist(similarity_scores, bins=num_bins)

    # Set colors!
    fracs = N / N.max()
    norm = colors.Normalize(fracs.min(), fracs.max())
    for thisfrac, thispatch in zip(fracs, patches):
        color = plt.cm.inferno(norm(thisfrac))
        thispatch.set_facecolor(color)

    # Set histogram title, axis labels, and ticks
    ax.set_title('Average similarity score for AMR genes across genome neighborhoods')
    ax.set_xlabel('Average similarity score')
    ax.set_ylabel('Number of neighborhoods')

    ytick_vals = [i + 1 for i in range(0, int(N.max()))]
    ax.set_yticks(ytick_vals)

    # Save figure to user specified output directory for clustering results
    savename = os.path.join(save_path + '/clustering/', 'avg_similarity_histogram.png')
    savefig(savename, bbox_inches='tight', dpi=100)


def plot_distance_histogram(distance_dict, save_path):
    """
    For use with clustering and scoring functionalities to plot a histogram of
    maximum differences measured across all AMR gene neighborhoods extracted.
    """
    distances = list(distance_dict.values())

    # Number of bins equal to the nearest multiple of 0.05 rounded up
    num_bins = np.arange(0.05, round(max(distances) / 0.05) * 0.05 + 0.05, 0.05)

    # Defining histogram
    fig, ax = plt.subplots(1, 1)
    N, bins, patches = ax.hist(distances, bins=num_bins)

    # Set colors
    fracs = N / N.max()
    norm = colors.Normalize(fracs.min(), fracs.max())
    for thisfrac, thispatch in zip(fracs, patches):
        color = plt.cm.inferno(norm(thisfrac))
        thispatch.set_facecolor(color)

    # Set histogram title, axis labels, and ticks
    ax.set_title('Average distance score obtained for AMR gene neighborhood across genomes')
    ax.set_xlabel('Average distance score')
    ax.set_ylabel('Number of neighborhoods')

    ax.set_xticks(num_bins)
    ytick_vals = [i + 1 for i in range(0, int(N.max()))]
    ax.set_yticks(ytick_vals)

    # Save figure to user specified output directory for clustering results
    savename = os.path.join(save_path + '/clustering/', 'mean_distances_histogram.png')
    savefig(savename, bbox_inches='tight', dpi=100)


def graph_UPGMA_clusters(linkage_matrix, labels, AMR_gene, output_path):
    """
    Visualizes UPGMA clusters using a dendrogram with branches coloured according to clusters.
    """
    fig = plt.figure(figsize=(15, 5))
    hierarchy.set_link_color_palette(['m', 'c', 'y', 'k'])
    hierarchy.dendrogram(Z=linkage_matrix, leaf_rotation=90, leaf_font_size=8, labels=labels)
    savename = os.path.join(output_path + '/clustering/UPGMA/', AMR_gene + ".png")
    plt.title("UPGMA dendrogram for gene {g} across {n} genomes".format(g=AMR_gene, n=len(labels)))
    plt.ylabel("Distance between neighborhoods")
    plt.xlabel("Neighborhoods with their respective Genome IDs")
    plt.savefig(savename, bbox_inches='tight', dpi=100)
    plt.close()


def plotly_dendrogram(linkage_matrix, labels, AMR_gene, output_path):
    """
    Generates an interactive dendrogram visualization using Plotly figure factory.
    """
    title = "UPGMA dendrogram for {g}".format(g=AMR_gene, n=len(labels))
    fig = create_dendrogram(linkage_matrix, labels=labels, colorscale=get_colors(len(labels)))
    fig.update_layout(autosize=True, title=title, paper_bgcolor='white', template='plotly_white', width=419, height=316)
    savename = os.path.join(output_path + '/clustering/UPGMA/', AMR_gene + ".html")
    fig.write_html(savename)
    # plotly.offline.plot(fig, filename=savename)


def draw_mcl_graph(matrix, clusters, **kwargs):
    """
    Modified version of draw_graph from Guy Allard to save image in specified directory rather than showing it.
    Original Source Code: https://github.com/GuyAllard/markov_clustering/blob/master/markov_clustering/drawing.py

    :param matrix: The unprocessed adjacency matrix
    :param clusters: list of tuples containing clusters as returned
                     by 'get_clusters'
    :param kwargs: Additional keyword arguments to be passed to
                   networkx.draw_networkx
    """
    # make a networkx graph from the adjacency matrix
    graph = nx.Graph(matrix)

    # map node to cluster id for colors
    cluster_map = {node: i for i, cluster in enumerate(clusters) for node in cluster}
    colors = [cluster_map[i] for i in range(len(graph.nodes()))]

    # if colormap not specified in kwargs, use a default
    if not kwargs.get("cmap", False):
        kwargs["cmap"] = cm.tab20

    # draw (beyond this point, original code was modified)
    fig = plt.figure()
    fig = nx.draw_networkx(graph, node_color=colors, **kwargs)
    axis("off")

    # save
    savename = os.path.join(kwargs['save_path'] + '/clustering/MCL/', kwargs['save_name'] + ".png")
    plt.savefig(savename, bbox_inches='tight', dpi=100)


def plotly_mcl_network(matrix, clusters, genome_names, AMR_gene, output_path):
    # make a networkx graph from the adjacency matrix
    graph = nx.Graph(matrix)
    pos = nx.spring_layout(graph)

    # Get node cluster assignments from MCL
    cluster_map = {node: i for i, cluster in enumerate(clusters) for node in cluster}

    # Determine edge lines between nodes: only add edges between nodes in the same cluster
    edge_x = []
    edge_y = []

    for cluster_group in clusters:
        cluster_nodes = list(cluster_group)
        for edge in graph.edges(cluster_nodes):
            if edge[0] in cluster_nodes and edge[1] in cluster_nodes:
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                edge_x.append(x0)
                edge_x.append(x1)
                edge_x.append(None)
                edge_y.append(y0)
                edge_y.append(y1)
                edge_y.append(None)

    # Assign a random color to every cluster
    hex_colors = []
    cluster_colors_dict = {}

    colors = get_colors(len(genome_names))
    for hex in colors:
        hex_colors.append(hex)

    for cluster in set(cluster_map.values()):
        cluster_colors_dict[str(cluster)] = hex_colors[cluster]

    # Get graph nodes
    node_x = []
    node_y = []
    for node in graph.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)

    # Add genome_ID as hover text to allow user to see which node comes from which genome neighborhood
    node_cluster = []
    node_text = {}

    node_to_genome_id_map = {}
    for node, genome_id in zip(cluster_map.keys(), genome_names):
        node_to_genome_id_map[node] = genome_id

    for i, cluster in enumerate(clusters):
        for node in cluster:
            node_cluster.append(i)
            node_text[node] = "{g}: Cluster {c}".format(g=node_to_genome_id_map[node], c=cluster_map[node])

    # Draw edges
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        showlegend=False,
        hoverinfo='none',
        mode='lines'
    )

    # Baseline figure with graph edges
    fig = go.Figure(data=edge_trace,
                    layout=go.Layout(
                        title='MCL network clusters for {}'.format(AMR_gene),
                        autosize=False,
                        width=419,
                        height=316,
                        showlegend=True,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )

    # Draw graph nodes: overlay on top of edge traces at respective positions
    legend_clusters = []
    for node in graph.nodes():
        x_pos, y_pos = pos[node]
        if cluster_map[node] not in legend_clusters:
            fig.add_trace(
                go.Scatter(x=[x_pos], y=[y_pos], name=str(cluster_map[node]), mode='markers',
                           hoverinfo='text', text=[node_text[node]], line=dict(width=2, color='black'),
                           legendgroup=str(cluster_map[node]), legendrank=cluster_map[node],
                           marker=dict(color=cluster_colors_dict[str(cluster_map[node])], size=15)))
            legend_clusters.append(cluster_map[node])

        else:
            fig.add_trace(
                go.Scatter(x=[x_pos], y=[y_pos], mode='markers',
                           hoverinfo='text', text=[node_text[node]],
                           legendgroup=str(cluster_map[node]), showlegend=False, legendrank=cluster_map[node],
                           marker=dict(color=cluster_colors_dict[str(cluster_map[node])], size=15)))

    fig.update_layout(legend_title='MCL Cluster', paper_bgcolor='white', template='plotly_white')

    savename = os.path.join(output_path + '/clustering/MCL/', AMR_gene + ".html")
    fig.write_html(savename)


def graph_DBSCAN_clusters(distance_matrix_df, DBSCAN_clusters, labels, AMR_gene, output_path):
    """
    Graphs DBSCAN clusters using Principal Coordinates Analysis (PCoA) plot.
    """
    # Get PCoA dimension values from distance matrix: preserve their relative relationships
    pcoa_vals = pcoa(distance_matrix_df)

    # Get number of unique clusters for color legend
    num_clusters = np.unique(labels)
    fig, ax = plt.subplots(1, 1)

    # Generalize color assignments to a generic number of points
    cmap = plt.cm.jet
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmap = cmap.from_list('custom cmap', cmaplist, cmap.N)

    bounds = np.linspace(0, len(num_clusters), len(num_clusters) + 1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    # Visualize PCoA as scatterplot with cluster labels as colors
    scatterplot = ax.scatter(x=pcoa_vals.samples['PC1'], y=pcoa_vals.samples['PC2'], c=labels, cmap=cmap, norm=norm)

    # Add genome labels to points
    for i in range(len(distance_matrix_df.index)):
        plt.text(pcoa_vals.samples.loc[str(i), 'PC1'], pcoa_vals.samples.loc[str(i), 'PC2'],
                 distance_matrix_df.index[i])

    # Create side colorbar legend
    cb = plt.colorbar(scatterplot, spacing='proportional', ticks=bounds)
    cb.set_label('DBSCAN cluster')

    # Title and axis labels
    plt.title("PCoA plot of {g} neighbourhoods".format(g=AMR_gene))
    plt.ylabel("PC1")
    plt.xlabel("PC2")

    savename = os.path.join(output_path + '/clustering/DBSCAN/', AMR_gene + ".png")
    plt.savefig(savename, bbox_inches='tight', dpi=100)


def plotly_pcoa(distance_matrix_df, genome_ids, labels, AMR_gene, output_path):
    """
    Make Plotly dash interactive scatterplot of PCoA visualization of DBSCAN clusters
    """
    # Get PCoA dimension values from distance matrix: preserve their relative relationships
    pcoa_vals = pcoa(distance_matrix_df)

    # Create new dataframe indexed by genome name containing data on PC1, PC2, and cluster label
    df_data = {'GenomeID': genome_ids,
               'PC1': pcoa_vals.samples['PC1'].values,
               'PC2': pcoa_vals.samples['PC2'].values,
               'Cluster': labels}

    df = pd.DataFrame(data=df_data, index=genome_ids)
    df['Cluster'] = df['Cluster'].astype(str)

    colors = get_colors(len(genome_ids))
    # Make noise cluster black by default
    colors.insert(0, '#111111')

    fig = px.scatter(df, x='PC1', y='PC2',
                     color='Cluster',
                     color_discrete_sequence=colors,
                     hover_name='GenomeID',
                     title='PCoA DBSCAN clusters for {g}'.format(g=AMR_gene))
    fig.update_traces(marker_size=5, line=dict(width=2, color='black'))
    fig.update_layout(paper_bgcolor='white', template='plotly_white', width=419, height=316)

    savename = os.path.join(output_path + '/clustering/DBSCAN/', AMR_gene + ".html")
    fig.write_html(savename)
