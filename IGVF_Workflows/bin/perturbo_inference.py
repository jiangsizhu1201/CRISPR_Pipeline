#!/usr/bin/env python

import argparse
import perturbo
import mudata as md
import numpy as np
import pandas as pd

def run_perturbo(mdata_input_fp, mdata_output_fp):
    """
    Run PerTurbo on the selected guide--element pairs and return a new MuData object
    with the test results stored in mdata.uns["test_results"]
    """
    mdata = md.read(mdata_input_fp)
    mdata["gene"].obs = (
        mdata.obs.join(mdata["gene"].obs.drop(columns=mdata.obs.columns, errors='ignore'))
        .join(mdata["guide"].obs.drop(columns=mdata.obs.columns.union(mdata["gene"].obs.columns), errors='ignore'))
        .assign(log1p_total_guide_umis=lambda x: np.log1p(x["total_guide_umis"]))
    )
    mdata["guide"].X = mdata["guide"].layers["guide_assignment"]
    pairs_to_test_df = pd.DataFrame(mdata.uns["pairs_to_test"])
    mdata.uns["intended_target_names"] = sorted(
        pd.unique(pairs_to_test_df["intended_target_name"])
    )
    # direct aggregate
    aggregated_df = (
        pairs_to_test_df.assign(value=1)
        .groupby(["gene_id", "intended_target_name"])
        .agg(value=('value', 'max')) 
        .reset_index()
    )

    # pivot the data
    mdata["gene"].varm["intended_targets"] = (
        aggregated_df
        .pivot(index="gene_id", columns="intended_target_name", values="value")
        .reindex(mdata["gene"].var_names)
        .fillna(0)
    )

    mdata.uns["intended_target_names"] = sorted(
        pd.unique(pairs_to_test_df["intended_target_name"])
    )

    intended_targets_df = pd.get_dummies(
        mdata["guide"].var["intended_target_name"]
    ).astype(float)

    mdata["guide"].varm["intended_targets"] = intended_targets_df[
        mdata.uns["intended_target_names"]
    ]
    perturbo.PERTURBO.setup_mudata(
            mdata,
            batch_key="batch",
            library_size_key="total_gene_umis",
            continuous_covariates_keys=["total_guide_umis"],
            guide_by_element_key="intended_targets",
            gene_by_element_key="intended_targets",
            modalities={
                "rna_layer": "gene",
                "perturbation_layer": "guide",
            },
        )

    model = perturbo.PERTURBO(mdata, likelihood="nb")
    model.train(20, lr=0.01, batch_size=128)

    igvf_name_map = {
        "element": "intended_target_name",
        "gene": "gene_id",
        "q_value": "p_value",
    }

    element_effects = (
        model.get_element_effects()
        .rename(columns=igvf_name_map)
        .assign(log2_fc=lambda x: x["loc"] / np.log(2))
        .merge(pairs_to_test_df)
    )

    mdata = md.read(mdata_input_fp)
    mdata.uns["test_results"] = element_effects[
        [
            "intended_target_name",
            "gene_id",
            "p_value",
            "log2_fc",
        ]
            ]
    
    mdata.write(mdata_output_fp)
    return mdata

def main():
    parser = argparse.ArgumentParser(description="Run PerTurbo analysis on MuData")
    parser.add_argument("mdata_input_fp", type=str, help="Input file path for MuData")
    parser.add_argument("mdata_output_fp", type=str, help="Output file path for MuData")

    args = parser.parse_args()
    run_perturbo(args.mdata_input_fp, args.mdata_output_fp)

if __name__ == "__main__":
    main()