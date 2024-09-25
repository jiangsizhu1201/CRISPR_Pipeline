#!/usr/bin/env python

import anndata as ad
import os
import re
import argparse
import pandas as pd

def extract_batch_num(filename):
    match = re.search(r'(.+)_ks_', filename)
    return match.group(1) if match else None

def process_adata_files(file_path_list, parsed_covariate_df):
    temp = pd.read_csv(parsed_covariate_df)
    cov_name = temp.columns.tolist()
    adatas = []

    # Sort the file path list by file name
    sorted_file_path_list = sorted(file_path_list, key=lambda x: os.path.basename(x))
    print(sorted_file_path_list)
    for file_path in sorted_file_path_list: 
        adata = ad.read_h5ad(os.path.join(file_path, "counts_unfiltered/adata.h5ad"))
        batch_num = extract_batch_num(os.path.basename(file_path))
        
        if batch_num:
            cov1 = cov_name[0]
            adata.obs[cov1] = batch_num
   
            adata.obs[cov1] = adata.obs[cov1].astype(str)
            temp[cov1] = temp[cov1].astype(str)
            
            adata.obs = adata.obs.join(temp.set_index(cov1), on=cov1)
        
        adatas.append(adata)
    
    return adatas

def main():
    parser = argparse.ArgumentParser(description='Process AnnData files and add batch information.')
    parser.add_argument('file_path_list', nargs='+', help='File path list for input AnnData files')
    parser.add_argument('parsed_covariate_df', type=str, help='Path to the parsed covariate DataFrame CSV')
    parser.add_argument('--output', type=str, default='combined_adata.h5ad', help='Output file name (default: combined_adata.h5ad)')

    args = parser.parse_args()

    processed_adatas = process_adata_files(args.file_path_list, args.parsed_covariate_df)
    combined_adata = ad.concat(processed_adatas, join='outer', index_unique="_")
    combined_adata.write_h5ad(args.output)
    print(f"Combined AnnData saved to {args.output}")

if __name__ == "__main__":
    main()
