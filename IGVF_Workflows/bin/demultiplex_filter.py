#!/usr/bin/env python

import anndata as ad
import pandas as pd
import argparse

def filter_adata_demux(adata_path, demux_report_path, demux_config_path, output_path):
    # Load the AnnData object
    adata = ad.read_h5ad(adata_path)
    demux_report = pd.read_csv(demux_report_path, index_col=0)
    demux_config = pd.read_csv(demux_config_path, header=None)
    demux_config.columns = ['cluster_id','hto_type']
    demux_config['hto_type'] = demux_config['hto_type'].str.strip()
    adata.obs['cluster_id'] = demux_report['Cluster_id']

    demux_config['hto_type_split'] = demux_config['hto_type'].str.split('-').str.join(",")
    demux_config['hto_type_split'] = demux_config['hto_type_split'].apply(lambda x: 'multiplets' if ',' in x else x)

    adata.obs = adata.obs.merge(
        demux_config,
        how='left',
        left_on='cluster_id',
        right_on='cluster_id'
        ).set_index(adata.obs.index)

    adata_filtered = adata[adata.obs['hto_type'] != 'negative'].copy()
    adata_filtered.write(output_path)

def main():
    # Set up argparse
    parser = argparse.ArgumentParser(description="Update AnnData object with demultiplexing report.")
    parser.add_argument('--adata', type=str, required=True, help="Path to the input h5ad file.")
    parser.add_argument('--demux_report', type=str, required=True, help="Path to the demultiplexing report CSV file.")
    parser.add_argument('--demux_config', type=str, required=True, help="Path to the demultiplexing report CSV file.")
    parser.add_argument('-o', '--output', type=str, required=True, help="Path to save the updated h5ad file.")

    args = parser.parse_args()
    filter_adata_demux(args.adata, args.demux_report, args.demux_config, args.output)

if __name__ == "__main__":
    main()
