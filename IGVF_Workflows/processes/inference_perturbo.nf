
process inference_perturbo {
    conda '/data/pinello/SHARED_SOFTWARE/anaconda_latest/envs/lb_envs/scverse_ml008'
    input:
    path input_mdata

    output:
    path "perturbo_mdata.h5mu", emit: inference_mudata

    script:
        """
        perturbo_inference.py ${input_mdata} perturbo_mdata.h5mu
        """
}
