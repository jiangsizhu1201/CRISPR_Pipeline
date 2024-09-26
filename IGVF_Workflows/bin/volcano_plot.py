#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import muon as mu
from mudata import MuData
import matplotlib.pyplot as plt
import os
from adjustText import adjust_text

def volcano(mdata: MuData, log2_fc_thresh: float, p_value_thresh: float):
    # Get test results from MuData
    test_results = pd.DataFrame({k: v for k, v in mdata.uns['test_results'].items()})
    
    # Plot all data points (non-significant genes)
    plt.scatter(x=test_results['log2_fc'], y=test_results['p_value'].apply(lambda x: -np.log10(x)), s=5, color="green", label="Not significant")
    
    # Highlight down or up-regulated genes based on the thresholds
    down = test_results[(test_results['log2_fc'] <= -log2_fc_thresh) & (test_results['p_value'] <= p_value_thresh)]
    up = test_results[(test_results['log2_fc'] >= log2_fc_thresh) & (test_results['p_value'] <= p_value_thresh)]

    filtered_down = down[np.isfinite(down['log2_fc']) & np.isfinite(down['p_value'])]

    print("Number of down-regulated genes based on the thresholds:", len(down))
    print("Number of down-regulated genes based on the thresholds (log2_fc and p_value are finite values):", len(filtered_down))
    print("Number of up-regulated genes based on the thresholds::", len(up))

    down_copy = down.copy()
    down_copy['log2_fc'] = down_copy['log2_fc'].replace([np.inf, -np.inf], [1e10, -1e10])
    down_copy['p_value'] = down_copy['p_value'].replace([np.inf, -np.inf], [1e10, -1e10])

    top_5_down_log2fc = down_copy.sort_values(by='log2_fc', ascending=True).head(10)
    top_5_down_pval = down_copy.sort_values(by='p_value', ascending=True).head(10)
    print("Annotated top 5 down-regulated genes based on the log2_fc thresholds:", top_5_down_log2fc)
    print("Annotated top 5 down-regulated genes based on the p value thresholds:", top_5_down_pval)
    annotated = pd.concat([top_5_down_log2fc, top_5_down_pval])

    # Plot down-regulated genes
    plt.scatter(x=down['log2_fc'], y=down['p_value'].apply(lambda x: -np.log10(x)), s=10, label="Down-regulated", color="blue")
    
    # Plot up-regulated genes
    plt.scatter(x=up['log2_fc'], y=up['p_value'].apply(lambda x: -np.log10(x)), s=10, label="Up-regulated", color="red")
    
    # Annotate top 10 down-regulated genes(sorted by log2_fc)
    texts = []
    # Annotate down-regulated genes 
    for index, row in annotated.iterrows():
        # Annotate with gene_id and guide_id, bold and size adjustment
        label = f"{row['intended_target_name']} ({row['guide_id']})"
        texts.append(plt.text(x=row['log2_fc'], y=-np.log10(row['p_value']), s=label, fontsize=8))

    # Adjust the text to prevent overlap
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='gray', lw=0.2))

    # Set plot boundaries and labels
    low, high = plt.xlim()
    bound = max(abs(low), abs(high)) + 0.5
    
    plt.xlabel("log2 fold change")
    plt.ylabel("p-value (log10)")
    
    # Draw threshold lines
    plt.axvline(-log2_fc_thresh, color="grey", linestyle="--")
    plt.axvline(log2_fc_thresh, color="grey", linestyle="--")
    plt.axhline(-np.log10(p_value_thresh), color="grey", linestyle="--")
    
    plt.xlim(-bound, bound)
    plt.legend(loc="upper right")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a volcano plot from MuData")
    parser.add_argument("mdata_path", type=str, help="Path to the MuData file")
    parser.add_argument("--log2_fc", type=float, default=1, help="log2 fold change threshold")
    parser.add_argument("--p_value", type=float, default=0.01, help="p-value threshold")
    
    args = parser.parse_args()
    
    # Load MuData
    mdata = mu.read(args.mdata_path)
    volcano(mdata, args.log2_fc, args.p_value)
    
    # Save the plot
    output_dir = "evaluation_output"
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "volcano_plot.png")

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
