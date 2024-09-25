#!/usr/bin/env python

import anndata as ad
import argparse
import os

def process_adata_files(file_paths):
    sorted_file_paths = sorted(file_paths, key=os.path.basename)
    print("Sorted file paths:", sorted_file_paths)  # Print sorted file paths
    return [ad.read_h5ad(path) for path in sorted_file_paths]

def main():
    parser = argparse.ArgumentParser(description='Process AnnData files and add batch information.')
    parser.add_argument('file_path_list', nargs='+', help='List of AnnData file paths')
    parser.add_argument('--output', type=str, required=True, help='Output file name')

    args = parser.parse_args()

    combined_adata = ad.concat(process_adata_files(args.file_path_list), join='outer')
    combined_adata.write_h5ad(args.output)
    print(f"Combined AnnData saved to {args.output}")

if __name__ == "__main__":
    main()
