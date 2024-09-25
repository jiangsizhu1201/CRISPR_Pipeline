#!/usr/bin/env python
import argparse
import asyncio
import muon as mu
from muon import MuData
import numpy as np
from scipy.sparse import dok_matrix
import cmdstanpy
import os

# Update these paths to point to compiled models
CS_MODEL_FILE = os.path.expanduser("~/IGVF_CRISPR_Pipeline/Workflows/IGVF_Workflows/cmdstan/cs-guide-mixture")
DC_MODEL_FILE = os.path.expanduser("~/IGVF_CRISPR_Pipeline/Workflows/IGVF_Workflows/cmdstan/dc-guide-mixture")

def posteriors_layer(stan_results, array, threshold=None):
    if threshold is None:
        for guide_id, (samples, cell_info) in stan_results.items():
            pzi = np.transpose(samples.stan_variable("PZi"))
            for i, (cell_id, *_) in enumerate(cell_info):
                array[cell_id, guide_id] = np.median(pzi[i])
    else:
        for guide_id, (samples, cell_info) in stan_results.items():
            pzi = np.transpose(samples.stan_variable("PZi"))
            for i, (cell_id, *_) in enumerate(cell_info):
                if np.median(pzi[i]) >= threshold:
                    array[cell_id, guide_id] = 1
    return array.tocsr()

async def run_cleanser(counts, model):
    stan_model = cmdstanpy.CmdStanModel(exe_file=model)
    
    N = len(counts)
    X = [count for _, _, count in counts]
    L = [1.0] * N  # Assuming L is a normalization factor, set to 1 for all observations
    
    stan_data = {
        "N": N,
        "X": X,
        "L": L
    }
    
    fit = stan_model.sample(
        data=stan_data, 
        iter_sampling=50,  # Reduced from 1000
        iter_warmup=50,    # Reduced from 1000
        chains=1,         
        parallel_chains=1  
    )
    
    results = {}
    for guide_id in set(guide for guide, _, _ in counts):
        cell_info = [(cell, count) for g, cell, count in counts if g == guide_id]
        results[guide_id] = (fit, cell_info)
    
    return results

def cleanser_posteriors(guides, threshold):
    guide_count_array = guides.X.todok()
    counts = [
        (key[1], key[0], int(guide_count))
        for key, guide_count in guide_count_array.items()
    ]
    analysis = guides.uns.get("capture_method")
    if analysis is None or analysis[0] == "CROP-seq":
        model = CS_MODEL_FILE
    elif analysis == "direct capture":
        model = DC_MODEL_FILE
    else:
        raise ValueError("Invalid capture method type")
    results = asyncio.run(run_cleanser(counts, model))
    return posteriors_layer(results, dok_matrix(guides.X.shape), threshold)

def threshold_posteriors(guides, threshold):
    guide_count_array = guides.X.todok()
    threshold = 5 if threshold is None else threshold
    array = dok_matrix(guides.X.shape)
    for (x, y), guide_count in guide_count_array.items():
        if guide_count >= threshold:
            array[x, y] = 1
    return array.tocsr()

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=str, help="Input MuData file")
    parser.add_argument("-t", "--threshold", default=None, type=float)
    model_group = parser.add_mutually_exclusive_group(required=True)
    model_group.add_argument(
        "--cleanser",
        action="store_true",
        help="Use CLEANSER to determine assignments")
    model_group.add_argument(
        "--umi-threshold",
        action="store_true",
        help="Use UMI threshold to determine assignments (default 5)")
    parser.add_argument("-o", "--output", required=True, type=str, help="Output MuData file")
    return parser.parse_args()

def run(gas, cleanser, umi_threshold, threshold):
    guides = gas["guide"]
    if cleanser:
        posteriors = cleanser_posteriors(guides, threshold)
    elif umi_threshold:
        posteriors = threshold_posteriors(guides, threshold)
    guides.layers["guide_assignment"] = posteriors
    guides.layers["guide_assignment"] = guides.layers["guide_assignment"].toarray()
    print(f"Shape of guide_assignment layer: {guides.layers['guide_assignment'].shape}")

if __name__ == "__main__":
    args = get_args()
    mu_input = mu.read(args.input)
    run(mu_input, args.cleanser, args.umi_threshold, args.threshold)
    mu.write(args.output, mu_input)