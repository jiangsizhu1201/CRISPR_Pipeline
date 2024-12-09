#!/usr/bin/env python

import argparse
import pandas as pd
import mudata as mu
import numpy as np

def export_output(mudata, inference_method):
    """
    Generate per-guide and per-element outputs from MuData and save as TSV files.
    """
    # Generate per-guide output with dynamic column renaming based on inference method
    col_map = {
        'sceptre': {'log2_fc': 'sceptre_log2_fc', 'p_value': 'sceptre_p_value'},
        'perturbo': {'log2_fc': 'perturbo_log2_fc', 'p_value': 'perturbo_p_value'}
    }
    if inference_method not in col_map:
        raise ValueError("Invalid inference_method. Must be 'sceptre' or 'perturbo'.")

    per_guide_output = (
        pd.DataFrame(mudata.uns['test_results'])
        .merge(mudata.mod['guide'].var, how='left', on='intended_target_name')
        .assign(
            cell_number=mudata.shape[0],
            avg_gene_expression=[
                np.mean(mudata.mod['gene'].X[:, mudata.mod['gene'].var['symbol'] == target].toarray()) 
                if (mudata.mod['gene'].var['symbol'] == target).any() 
                else np.nan
                for target in pd.DataFrame(mudata.uns['test_results'])['intended_target_name']
            ]
        )
        .rename(columns=col_map[inference_method])
    )

    # Ensure missing columns are added
    for col in ['sceptre_log2_fc', 'sceptre_p_value', 'perturbo_log2_fc', 'perturbo_p_value']:
        if col not in per_guide_output:
            per_guide_output[col] = None

    # Reorder columns
    per_guide_output = per_guide_output[
        ['intended_target_name', 'guide_id_x', 'intended_target_chr', 
         'intended_target_start', 'intended_target_end', 'gene_id', 
         'sceptre_log2_fc', 'sceptre_p_value', 'perturbo_log2_fc', 
         'perturbo_p_value', 'cell_number', 'avg_gene_expression']
    ].rename(columns={'guide_id_x': 'guide_id(s)'})

    # Generate per-element output
    per_element_output = (
        per_guide_output.groupby('intended_target_name', as_index=False)
        .agg({'guide_id(s)': ','.join})
        .merge(
            per_guide_output.drop_duplicates('intended_target_name')
            .drop(columns=['guide_id(s)']),
            on='intended_target_name',
            how='left'
        )
    )

    # Save outputs
    per_guide_output.to_csv('per_guide_output.tsv', sep='\t', index=False)
    per_element_output.to_csv('per_element_output.tsv', sep='\t', index=False)
    print("Exported output files successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Export outputs from MuData.")
    parser.add_argument("--inference_method", type=str, required=True, help="Inference method: 'sceptre' or 'perturbo'.")
    parser.add_argument("--mudata", type=str, required=True, help="Path to input MuData file.")
    args = parser.parse_args()

    # Load MuData and call the export function
    mudata = mu.read_h5mu(args.mudata)
    export_output(mudata, args.inference_method)
