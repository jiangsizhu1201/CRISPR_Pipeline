#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from GTFProcessing import GTFProcessing
import muon as mu

def main(limit, mudata_path, input_gtf_path):
    mudata = mu.read(mudata_path)
    guide_data = mudata.mod['guide'].var
    
    # Load GTF
    gtf = GTFProcessing(input_gtf_path)
    df_gtf_refseq = gtf.get_gtf_df()
    df_gtf_unique_gene_name = df_gtf_refseq.drop_duplicates('gene_name')
    df_gtf_unique_gene_name.set_index('gene_name', inplace=True)

    final_candidates_list = []

    for g_name, target_name, g_chr, g_start in guide_data[['guide_id', 'intended_target_name', 'intended_target_chr', 'intended_target_start']].values:
        temp_gtf = df_gtf_unique_gene_name.query(f'chr == "{g_chr}"').copy()
        temp_gtf['distance_from_guide'] = np.abs(temp_gtf['start'] - g_start)
        final_candidates = temp_gtf.query(f'distance_from_guide <= {limit}').copy()
        final_candidates['guide_id'] = g_name
        final_candidates['intended_target_name'] = target_name
        final_candidates['pair_type'] = 'discovery'
        final_candidates.reset_index(inplace=True)
        final_candidates_list.append(final_candidates[['guide_id', 'gene_name', 'intended_target_name', 'pair_type']].copy())
    final_candidates_df = pd.concat(final_candidates_list, ignore_index=True)
    final_candidates_df.to_csv('pairs_to_test.csv', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process GTF and guide data.")
    parser.add_argument('--limit', type=int, default=1000000, help='Limit for distance from guide')
    parser.add_argument('mudata', type=str, help='Path to the input MuData file')
    parser.add_argument('input_gtf', type=str, help='Path to the input GTF file')

    args = parser.parse_args()
    main(args.limit, args.mudata, args.input_gtf)
