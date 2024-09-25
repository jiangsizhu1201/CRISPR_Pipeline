#!/usr/bin/env python

import argparse
import anndata as ad
import pandas as pd
import numpy as np
from muon import MuData
import matplotlib.pyplot as plt
import os
from GTFProcessing import GTFProcessing

def main(adata_rna, adata_guide, guide_metadata, gtf, moi):
    # Load the data
    guide_metadata = pd.read_excel(guide_metadata)
    adata_rna = ad.read_h5ad(adata_rna)
    adata_guide = ad.read_h5ad(adata_guide)
    gtf = GTFProcessing(gtf)
    df_gtf = gtf.get_gtf_df()   

    ## change in adata_guide
    # adding var for guide
    adata_guide.var.reset_index(inplace=True)
    # rename gene_id in guide
    adata_guide.var.rename(columns={'gene_id': 'guide_id'}, inplace=True)

     # check if the lengths are the same
    if len(guide_metadata) != len(adata_guide.var):
        print(f"The numbers of sgRNA_ID/guide_id are different: There are {len(adata_guide.var)} in guide anndata, but there are {len(guide_metadata)} in guide metadata.")

    # guide_metadata['guide_number'] = guide_metadata['sgRNA_ID'].apply(lambda x: 1 if 'sg1' in x else (2 if 'sg2' in x else None))
    adata_guide.var = adata_guide.var.merge(guide_metadata, left_on='guide_id', right_on ='sgRNA_ID', how='inner')
    adata_guide.var[['intended_target_name', 'intended_target_chr', 'intended_target_start', 'intended_target_end', 'sequence']] = guide_metadata[['Target_name', 'chr', 'start', 'end', 'sgRNA_sequences']]
    adata_guide.var['targeting'] =  'TRUE'

    # reset feature_id to var_names (index)
    guide_metadata['feature_id'] = guide_metadata['sgRNA_ID'] + "|" + guide_metadata['sgRNA_sequences']
    adata_guide.var_names = guide_metadata['feature_id']

    # adding uns for guide (needed for sceptre_inference; moi as low/high)
    adata_guide.uns['moi'] = np.array([moi], dtype=object)
    adata_guide.uns['capture_method'] = np.array(['CROP-seq'], dtype=object)

    # adding number of nonzero guides and batch number
    adata_guide.obs['num_expressed_guides'] = (adata_guide.X > 0).sum(axis=1)
    adata_guide.obs['batch_number'] = adata_guide.obs['batch'].factorize()[0] + 1
    adata_guide.obs['total_guide_umis'] = adata_guide.X.sum(axis=1)
    
    ## change in adata_rna; 
    df_gtf['gene_id2'] = df_gtf['gene_id'].str.split('.').str[0]
    df_gtf = df_gtf.drop_duplicates('gene_id2')
    df_gtf_copy = df_gtf.copy()
    df_gtf_copy.set_index('gene_id2', inplace=True)
    # adding gene_start, gene_end, gene_chr
    adata_rna.var = adata_rna.var.join(df_gtf_copy[['chr', 'start', 'end']].rename(columns={'chr': 'gene_chr', 
                                                                                            'start': 'gene_start', 
                                                                                            'end': 'gene_end'}))
 
    # rename adata_rna obs
    adata_rna.obs.rename(columns={'n_genes_by_counts': 'n_counts',
                                'pct_counts_mt': 'percent_mito',
                                'n_genes' : 'num_expressed_genes',
                                'total_counts' : 'total_gene_umis'}, inplace=True)
    

    # knee plots
    knee_df = pd.DataFrame({
        'sum': np.array(adata_guide.X.sum(1)).flatten(),
        'barcodes': adata_guide.obs_names.values})
    knee_df = knee_df.sort_values('sum', ascending=False).reset_index(drop=True)
    knee_df['sum_log'] = np.log1p(knee_df['sum'])

    plt.figure(figsize=(8, 5))
    plt.plot(knee_df.index, knee_df['sum_log'], marker='o', linestyle='-', markersize=3)
    plt.xlabel('Barcode Index')
    plt.ylabel('Log of UMI Counts')
    plt.title('Knee Plot')

    # Save knee plot
    if not os.path.exists('figures'):
        os.makedirs('figures')
        print(f"Directory '{'figures'}' created.")
    else:
        print(f"Directory already exists.")

    plt.savefig('figures/knee_plot_guide.png')

    # Find the intersection of barcodes between scRNA and guide data
    intersecting_barcodes = list(set(adata_rna.obs_names)
                                 .intersection(adata_guide.obs_names))

    mdata = MuData({
        'gene': adata_rna[intersecting_barcodes, :].copy(),
        'guide': adata_guide[intersecting_barcodes, :].copy()
    })

    obs_names = set(mdata.mod['guide'].obs.columns.tolist()) & set(mdata.mod['gene'].obs.columns.tolist())
    mdata.obs = mdata.mod['guide'].obs.loc[:, list(obs_names)]

    # Save the MuData object
    mdata.write("mudata.h5mu")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a MuData object from scRNA and guide data.')
    parser.add_argument('adata_rna', type=str, help='Path to the scRNA AnnData file.')
    parser.add_argument('adata_guide', type=str, help='Path to the guide AnnData file.')
    parser.add_argument('guide_metadata', type=str, help='Path to the guide metadata Excel file.')
    parser.add_argument('gtf', type=str, help='Path to the GTF file.')
    parser.add_argument('moi', type=str, help='MOI information.')
    

    args = parser.parse_args()
    main(args.adata_rna, args.adata_guide, args.guide_metadata, args.gtf, args.moi)
