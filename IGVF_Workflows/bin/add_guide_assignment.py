#!/usr/bin/env python
import muon as mu
import pandas as pd
import numpy as np
import argparse
from scipy.sparse import csr_matrix

def add_guide_assignment(mudata, guide_assignment_csv):
    # Load MuData object
    mudata = mu.read_h5mu(mudata)

    # Read and process guide assignment
    df = pd.read_csv(guide_assignment_csv)
    sparse_matrix = csr_matrix(df.T)
    mudata.mod['guide'].layers['guide_assignment'] = sparse_matrix

    # Save the modified MuData object
    mudata.write("sceptre_assignment_mudata.h5mu")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process MuData and guide assignment.')
    parser.add_argument('--mudata', required=True, help='Path to the input h5mu file.')
    parser.add_argument('--guide_assignment_csv', required=True, help='Path to the guide assignment CSV file.')
    
    args = parser.parse_args()
    add_guide_assignment(args.mudata, args.guide_assignment_csv)