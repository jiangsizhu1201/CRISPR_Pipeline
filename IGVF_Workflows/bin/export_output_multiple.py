#!/usr/bin/env python

import argparse
import pandas as pd
import mudata as mu
import numpy as np 

def merge_results(sceptre_result, perturbo_mudata):

    # Load and rename SCEPTRE results
    test_result = pd.read_csv(sceptre_result).rename(
        columns={'log2_fc': 'sceptre_log2_fc', 'p_value': 'sceptre_p_value'}
    )

    # Load and rename Perturbo results
    mudata = mu.read_h5mu(perturbo_mudata)
    perturbo_result = pd.DataFrame(mudata.uns['test_results']).rename(
        columns={'log2_fc': 'perturbo_log2_fc', 'p_value': 'perturbo_p_value'}
    )

    # Merge results 
    merged_result = test_result.merge(
        perturbo_result[['guide_id', 'perturbo_log2_fc', 'perturbo_p_value']],
        on='guide_id', how='outer'
    )

    # Update MuData with merged results
    mudata.uns["test_results"] = merged_result.to_dict(orient='list')
    mudata.write("inference_mudata.h5mu")
    print("Updated MuData saved to 'inference_mudata.h5mu'.")

    return merged_result, mudata

def export_output(merged_result, mudata):
    """
    Generate and save per-guide and per-element outputs.
    """
    # Generate per-guide output
    per_guide_output = (
        merged_result
        .merge(mudata.mod['guide'].var, how='left', on='intended_target_name')
        .assign(
            cell_number=mudata.shape[0],
            avg_gene_expression=[
                np.mean(mudata.mod['gene'].X[:, mudata.mod['gene'].var['symbol'] == target].toarray()) 
                if (mudata.mod['gene'].var['symbol'] == target).any() 
                else np.nan
                for target in merged_result['intended_target_name']
            ]
        )
        .rename(columns={'guide_id_x': 'guide_id(s)'})
        .reindex(columns=[
            'intended_target_name', 'guide_id(s)', 'intended_target_chr',
            'intended_target_start', 'intended_target_end', 'gene_id',
            'sceptre_log2_fc', 'sceptre_p_value', 'perturbo_log2_fc',
            'perturbo_p_value', 'cell_number', 'avg_gene_expression'
        ])
    )
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
    parser = argparse.ArgumentParser(description='Process and merge SCEPTRE and Perturbo results.')
    parser.add_argument('--sceptre_result', required=True, help='Path to the SCEPTRE results CSV file.')
    parser.add_argument('--perturbo_mudata', required=True, help='Path to the Perturbo h5mu file.')
    args = parser.parse_args()

    # Process results and export outputs
    merged_result, mudata = merge_results(args.sceptre_result, args.perturbo_mudata)
    export_output(merged_result, mudata)
