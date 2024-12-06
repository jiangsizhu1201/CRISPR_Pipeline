#!/usr/bin/env python
import os
import argparse
import pandas as pd
import numpy as np
import muon as mu
from mudata import MuData
import networkx as nx
import matplotlib.pyplot as plt
from typing import Optional, Literal

def plot_network(mdata: MuData, central_node: str, method: Literal['sceptre', 'perturbo'],
                source_column: str = "source", target_column: str = "target",
                min_weight: Optional[float] = None, results_key: Optional[str] = "test_results",
                ax: Optional[plt.Axes] = None):
    
    if ax is None:
        ax = plt.gca()
    
    # Set method-specific column names
    weight_column = f"{method}_log2_fc"
    node_size_column = f"{method}_p_value"
    
    results_df = pd.DataFrame({k: v for k, v in mdata.uns[results_key].items()})
    results_df = results_df.drop_duplicates()
    if min_weight is not None:
        results_df = results_df[results_df[weight_column].abs() >= min_weight]

    # Filter to rows related to the selected central node
    if central_node is not None:
        results_df = results_df[results_df[source_column] == central_node]

    G = nx.DiGraph()
    for i, row in results_df.iterrows():
        G.add_edge(row[source_column], row[target_column],
                  weight=row[weight_column],
                  pvalue=row[node_size_column])
    pos = nx.circular_layout(G)

    # Calculate node sizes based on -log10(p-value)
    node_sizes = {}
    for node in G.nodes():
        edges = G.in_edges(node, data=True)
        if edges:
            min_pvalue = min(d['pvalue'] for _, _, d in edges)
            node_sizes[node] = -np.log10(min_pvalue) * 100
        else:
            node_sizes[node] = 100

    # Draw nodes
    sizes = [node_sizes.get(n, 100) for n in G.nodes()]
    nx.draw(G, pos, with_labels=True, node_color="skyblue",
            node_size=sizes, edge_cmap=plt.cm.coolwarm, arrowsize=20)

    label_pos = {k: (v[0]+0.2, v[1] + 0.05) for k, v in pos.items()}
    nx.draw_networkx_labels(G, label_pos, font_size=8)
    
    edge_labels = {(u, v): f"FC:{d['weight']:.2f}\np:{d['pvalue']:.2e}"
                  for u, v, d in G.edges(data=True)}
    edge_colors = [d["weight"] for u, v, d in G.edges(data=True)]
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors,
                          edge_cmap=plt.cm.coolwarm, arrowsize=5)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=6)
    
    ax.set_title(f"{method.capitalize()} Network - {central_node}")

def select_central_nodes(mdata: MuData, num_nodes: int, source_column: str,
                        target_column: str, method: Literal['sceptre', 'perturbo']):
    weight_column = f"{method}_log2_fc"
    
    results_df = pd.DataFrame({k: v for k, v in mdata.uns["test_results"].items()})
    results_df = results_df.drop_duplicates()
    
    # Create a graph and calculate the degree of each node
    G = nx.DiGraph()
    for i, row in results_df.iterrows():
        if pd.notna(row.get(weight_column)):  # Only add edges if weight exists
            G.add_edge(row[source_column], row[target_column], weight=abs(row[weight_column]))
    
    # Sort nodes by weighted degree
    degrees = dict(G.degree(weight='weight'))
    sorted_nodes = sorted(degrees, key=degrees.get, reverse=True)

    intended_target_names = set(results_df['intended_target_name'])
    filtered_nodes = [node for node in sorted_nodes if node in intended_target_names]
    print(f"Top {num_nodes} nodes for {method}:", filtered_nodes[:num_nodes])
    return filtered_nodes[:num_nodes]

def main():
    parser = argparse.ArgumentParser(description="Select nodes and plot networks comparing two analysis methods")
    parser.add_argument("mdata_path", type=str, help="Path to the MuData file")
    parser.add_argument("--num_nodes", type=int, required=True, help="Number of central nodes to select")
    parser.add_argument("--source_column", type=str, default="intended_target_name",
                      help="Source column name in test results")
    parser.add_argument("--target_column", type=str, default="gene_id",
                      help="Target column name in test results")
    parser.add_argument("--min_weight", type=float, default=0.1,
                      help="Minimum absolute log2 fold change to filter edges")
    parser.add_argument("--results_key", type=str, default="test_results",
                      help="Key for test results in mdata.uns")
    
    args = parser.parse_args()
    
    # Load MuData
    mdata = mu.read(args.mdata_path)

    # Select central nodes for both methods
    sceptre_nodes = set(select_central_nodes(
        mdata, args.num_nodes, args.source_column, args.target_column, 'sceptre'
    ))
    perturbo_nodes = set(select_central_nodes(
        mdata, args.num_nodes, args.source_column, args.target_column, 'perturbo'
    ))
    
    # Union of top nodes from both methods
    central_nodes = list(sceptre_nodes | perturbo_nodes)
    print(f"\nTotal unique nodes selected: {len(central_nodes)}")

    # Create a grid of plots
    fig = plt.figure(figsize=(15, 7.5 * len(central_nodes)))
    
    for i, central_node in enumerate(central_nodes):
        # Plot Sceptre
        ax1 = fig.add_subplot(len(central_nodes), 2, 2*i + 1)
        plot_network(
            mdata,
            central_node=central_node,
            method='sceptre',
            source_column=args.source_column,
            target_column=args.target_column,
            min_weight=args.min_weight,
            results_key=args.results_key,
            ax=ax1
        )
        
        # Plot Perturbo
        ax2 = fig.add_subplot(len(central_nodes), 2, 2*i + 2)
        plot_network(
            mdata,
            central_node=central_node,
            method='perturbo',
            source_column=args.source_column,
            target_column=args.target_column,
            min_weight=args.min_weight,
            results_key=args.results_key,
            ax=ax2
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
