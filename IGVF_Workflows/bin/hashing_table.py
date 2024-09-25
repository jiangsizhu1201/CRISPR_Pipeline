#!/usr/bin/env python

import pandas as pd
import argparse

def process_table(hashing_table):
    hashing_metadata = pd.read_excel(hashing_table, header=0)
    hashing_metadata[['HTO_sequences', 'HTO_ID']].to_csv('hashing_table.txt',
                                                sep='\t', header=None, index=None)

def main():
    parser = argparse.ArgumentParser(description="Process hashing metadata.")
    parser.add_argument('--hashing_table', type=str, required=True, help="Path to the input hashing file.")
    args = parser.parse_args()

    process_table(args.hashing_table)

if __name__ == "__main__":
    main()
