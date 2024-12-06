#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
import muon as mu
import os
from GTFProcessing import GTFProcessing
from typing import Literal, Tuple, Dict

def process_coordinates(mdata) -> Dict[str, list]:
    """Extract coordinate information from MuData object"""
    coordinate_dict = {}
    
    # Process gene coordinates
    for index, row in mdata.mod["gene"].var.iterrows():
        if np.isnan(row["gene_start"]) or np.isnan(row["gene_end"]):
            continue
        coordinate_dict[index] = [row["gene_chr"], row["gene_start"], row["gene_end"]]

    # Process guide coordinates
    for index, row in mdata.mod["guide"].var.iterrows():
        if row["intended_target_name"] in coordinate_dict or row["intended_target_name"] == "non-targeting":
            continue
        coordinate_dict[row["intended_target_name"]] = [
            row["intended_target_chr"],
            row["intended_target_start"],
            row["intended_target_end"]
        ]
    
    return coordinate_dict

def igv(mdata, gtf: str, method: Literal['sceptre', 'perturbo']) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Generate bedpe and bedgraph data for a specific method"""
    # Set method-specific column names
    log2_fc_col = f"{method}_log2_fc"
    p_value_col = f"{method}_p_value"
    
    # Process coordinates
    coordinate_dict = process_coordinates(mdata)

    # Process GTF file
    gtf_processor = GTFProcessing(gtf)
    df_gtf = gtf_processor.get_gtf_df()
    gencode_df = df_gtf[['gene_id', 'gene_name']]
    gencode_df['gene_id2'] = gencode_df['gene_id'].str.split('.').str[0]
    gencode_df = gencode_df.drop_duplicates()

    # Initialize data structures
    bedpe = defaultdict(list)
    bedgraph = defaultdict(list)
    
    # Process test results
    test_results = pd.DataFrame({k: v for k, v in mdata.uns['test_results'].items()})
    merged_df = test_results.merge(
        gencode_df[['gene_id2', 'gene_name']], 
        left_on='gene_id', 
        right_on='gene_id2', 
        how='left'
    )

    # Filter out rows where method-specific columns are missing
    merged_df = merged_df.dropna(subset=[log2_fc_col, p_value_col])

    for index, row in merged_df.iterrows():
        if row["intended_target_name"] == row["gene_name"]:
            # PROMOTER interactions
            if row["intended_target_name"] in coordinate_dict:
                coords = coordinate_dict[row["intended_target_name"]]
                bedgraph["chr"].append(coords[0])
                bedgraph["start"].append(coords[1])
                bedgraph["end"].append(coords[2])
                bedgraph["p_value"].append(row[p_value_col])
                bedgraph["log2_fc"].append(row[log2_fc_col])
        else:
            # ENHANCER-GENE interactions
            if (row["intended_target_name"] in coordinate_dict and 
                row["gene_id"] in coordinate_dict):
                source_coords = coordinate_dict[row["intended_target_name"]]
                target_coords = coordinate_dict[row["gene_id"]]
                
                bedpe["chr1"].append(source_coords[0])
                bedpe["start1"].append(source_coords[1])
                bedpe["end1"].append(source_coords[2])
                bedpe["chr2"].append(target_coords[0])
                bedpe["start2"].append(target_coords[1])
                bedpe["end2"].append(target_coords[2])
                bedpe["p_value"].append(row[p_value_col])
                bedpe["log2_fc"].append(row[log2_fc_col])

    bedpe_df = pd.DataFrame(bedpe)
    bedgraph_df = pd.DataFrame(bedgraph)

    if bedpe_df.empty:
        print(f"Warning: {method} bedpe_df is empty.")
    if bedgraph_df.empty:
        print(f"Warning: {method} bedgraph_df is empty.")
    
    print(f"\n{method.capitalize()} statistics:")
    print(f"Number of enhancer-gene interactions: {len(bedpe_df)}")
    print(f"Number of promoter interactions: {len(bedgraph_df)}")

    return bedpe_df, bedgraph_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process MuData and generate bedpe and bedgraph files for multiple methods")
    parser.add_argument("mdata_path", type=str, help="Path to the MuData file")
    parser.add_argument("--gtf", type=str, required=True, help="Path to the GTF file")
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = "evaluation_output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load MuData
    mdata = mu.read(args.mdata_path)
    
    # Process data for both methods
    for method in ['sceptre', 'perturbo']:
        bedpe_df, bedgraph_df = igv(mdata, args.gtf, method)
        
        # Generate outputs with fixed filenames
        bedpe_path = os.path.join(output_dir, f"{method}.bedpe")
        bedgraph_path = os.path.join(output_dir, f"{method}.bedgraph")
        
        # Save files
        bedpe_df.to_csv(bedpe_path, sep="\t", index=False, header=False)
        bedgraph_df.to_csv(bedgraph_path, sep="\t", index=False, header=False)
        
        print(f"\n{method.capitalize()} files saved:")
        print(f"bedpe file: {bedpe_path}")
        print(f"bedgraph file: {bedgraph_path}")
