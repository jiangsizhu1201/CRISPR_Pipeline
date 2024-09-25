#!/usr/bin/env python
import os
import argparse
import pandas as pd
import numpy as np
import muon as mu
from mudata import MuData

import networkx as nx
import matplotlib.pyplot as plt
from typing import Optional

# Function to plot a network
def plot_network(mdata: MuData, central_node: str, source_column: str = "source",
                 target_column: str = "target", weight_column: Optional[str] = None, min_weight: Optional[float] = None,
                 node_size_column: Optional[str] = None, results_key: Optional[str] = "test_results", ax: Optional[plt.Axes] = None):
    
    if ax is None:
        ax = plt.gca()
    
    results_df = pd.DataFrame({k: v for k, v in mdata.uns[results_key].items()})
    results_df = results_df.drop_duplicates()
    if min_weight is not None:
        results_df = results_df[results_df[weight_column].abs() >= min_weight]

    # Filter to rows related to the selected central node
    if central_node is not None:
        results_df = results_df[results_df[source_column] == central_node]

    G = nx.DiGraph()
    for i, row in results_df.iterrows():
        G.add_edge(row[source_column], row[target_column], weight=row[weight_column])
    pos = nx.circular_layout(G)

    if node_size_column is not None:
        node_size = results_df.set_index(target_column)[node_size_column].to_dict()
        sizes = [node_size.get(n, 10) for n in G.nodes]
        sizes = [s * 100 for s in sizes]
        nx.draw(G, pos, with_labels=True, node_color="skyblue", node_size=sizes, edge_cmap=plt.cm.Blues, arrowsize=20)
    else:
        nx.draw(G, pos, with_labels=True, node_color="skyblue", node_size=150, edge_cmap=plt.cm.Blues, arrowsize=20)

    label_pos = {k: (v[0]+0.2, v[1] + 0.05) for k, v in pos.items()}
    nx.draw_networkx_labels(G, label_pos, font_size=8)
    
    edge_labels = {(u, v): f"{d['weight']:.2f}" for u, v, d in G.edges(data=True)}
    nx.draw_networkx_edges(G, pos, edge_color=[d["weight"] for u, v, d in G.edges(data=True)], edge_cmap=plt.cm.coolwarm, arrowsize=5)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)

# Function to select central nodes based on degree
def select_central_nodes(mdata: MuData, num_nodes: int, source_column: str, target_column: str, weight_column: Optional[str] = None):
    results_df = pd.DataFrame({k: v for k, v in mdata.uns["test_results"].items()})
    results_df = results_df.drop_duplicates()
    # Create a graph and calculate the degree of each node
    G = nx.DiGraph()
    for i, row in results_df.iterrows():
        G.add_edge(row[source_column], row[target_column], weight=row.get(weight_column, 1))
    
    # Sort nodes by degree and select the top num_nodes
    degrees = dict(G.degree(weight=weight_column))
    sorted_nodes = sorted(degrees, key=degrees.get, reverse=True)

    intended_target_names = set(results_df['intended_target_name'])
    filtered_nodes = [node for node in sorted_nodes if node in intended_target_names]
    print(filtered_nodes[:num_nodes])
    return filtered_nodes[:num_nodes]

# Main function
def main():
    parser = argparse.ArgumentParser(description="Select nodes and plot a network from MuData")
    parser.add_argument("mdata_path", type=str, help="Path to the MuData file")
    parser.add_argument("--num_nodes", type=int, required=True, help="Number of central nodes to select")
    parser.add_argument("--source_column", type=str, default="intended_target_name", help="Source column name in test results")
    parser.add_argument("--target_column", type=str, default="gene_id", help="Target column name in test results")
    parser.add_argument("--weight_column", type=str, default="log2_fc", help="Weight column name in test results")
    parser.add_argument("--min_weight", type=float, default=0.1, help="Minimum weight to filter edges")
    parser.add_argument("--node_size_column", type=str, default="p_value", help="Column to determine node size")
    parser.add_argument("--results_key", type=str, default="test_results", help="Key for test results in mdata.uns")
    
    args = parser.parse_args()
    
    # Load MuData
    mdata = mu.read(args.mdata_path)

    # Select central nodes
    central_nodes = select_central_nodes(
        mdata,
        num_nodes=args.num_nodes,
        source_column=args.source_column,
        target_column=args.target_column,
        weight_column=args.weight_column
    )

    # Create a grid of plots
    cols = 2  
    rows = (len(central_nodes) + cols - 1) // cols

    fig = plt.figure(figsize=(15, 7.5 * rows))
    
    for i, central_node in enumerate(central_nodes):
        ax = fig.add_subplot(rows, cols, i + 1)
        plot_network(
            mdata,
            central_node=central_node,
            source_column=args.source_column,
            target_column=args.target_column,
            weight_column=args.weight_column,
            min_weight=args.min_weight,
            node_size_column=args.node_size_column,
            results_key=args.results_key,
            ax=ax
        )
    plt.tight_layout()

    # Save the plot
    output_dir = "evaluation_output"
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "network_plot.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_file}")

if __name__ == "__main__":
    main()
