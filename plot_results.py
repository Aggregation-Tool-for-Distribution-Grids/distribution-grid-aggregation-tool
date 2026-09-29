import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx

import numpy as np


def plot_results(
        orig_adj,
        final_adj,
        num_nodes,
        slack_buses,
        active_supernodes,
        clusters
):

    ### plot
    fig, (ax1, ax2) = plt.subplots(
        1,
        2,
        figsize=(16, 7)
    )

    ## original network graph
    G_orig = nx.Graph()
    rows_orig, cols_orig = np.where(orig_adj == 1)

    for i in range(num_nodes):
        G_orig.add_node(i)

    for u, v in zip(rows_orig, cols_orig):
        if u < v:
            G_orig.add_edge(u, v)

    # remove nodes that not connected
    G_orig.remove_nodes_from(list(nx.isolates(G_orig)))

    # determining node coordinates
    pos = nx.spring_layout(G_orig, seed=42)

    # colours for nodes
    # - red: slack node
    # - blue: others

    orig_node_colors = []
    for node in G_orig.nodes():
        if node in slack_buses:
            orig_node_colors.append('red')

        else:
            orig_node_colors.append('blue')

    # draw original graph

    nx.draw_networkx_nodes(
        G_orig,
        pos,
        ax=ax1,
        node_color=orig_node_colors,
        node_size=150,
        alpha=0.25
    )

    nx.draw_networkx_edges(
        G_orig,
        pos,
        ax=ax1,
        width=1.5,
        edge_color='gray'
    )

    nx.draw_networkx_labels(
        G_orig,
        pos,
        ax=ax1,
        font_size=5
    )

    ## graph after Opti-KRON
    G_red = nx.Graph()

    for sn in active_supernodes:
        G_red.add_node(sn)

    red_node_colors = []
    rows_red, cols_red = np.where(final_adj == 1)

    for u, v in zip(rows_red, cols_red):
        if u < v and u in active_supernodes and v in active_supernodes:
            G_red.add_edge(u, v)

    # remove nodes that not connected
    G_red.remove_nodes_from(list(nx.isolates(G_red)))

    # determining node coordinates
    pos = nx.spring_layout(G_red)

    red_node_colors = []
    for node in G_red.nodes():
        if node in slack_buses:
            red_node_colors.append('red')

        elif len(clusters[node]) > 1:
            # super node
            red_node_colors.append('orange')

        else:
            red_node_colors.append('blue')

    nx.draw_networkx_nodes(
        G_red,
        pos,
        ax=ax2,
        nodelist=list(G_red.nodes()),
        node_color=red_node_colors,
        node_size=150,
        alpha=0.25
    )

    nx.draw_networkx_edges(
        G_red,
        pos,
        ax=ax2,
        width=1.5,
        edge_color='gray'
    )

    nx.draw_networkx_labels(
        G_red,
        pos,
        ax=ax2,
        font_size=5
    )

    ax1.set_title("Original Network")
    ax2.set_title("After Opti-KRON")

    legend_elements = [
        Line2D(
            [0],
            [0],
            marker='o',
            color='w',
            label='Slack Bus',
            markerfacecolor='red',
            markersize=10
        ),
        Line2D(
            [0],
            [0],
            marker='o',
            color='w',
            label='Supernode (Merged)',
            markerfacecolor='orange',
            markersize=10
        ),
        Line2D(
            [0],
            [0],
            marker='o',
            color='w',
            label='Standard Bus',
            markerfacecolor='blue',
            markersize=10
        )
    ]

    ax2.legend(
        handles=legend_elements,
        loc='lower center',
        frameon=True,
        bbox_to_anchor=(-0.1, -0.1),
        ncol=3
    )

    plt.show()