
process inference_perturbo {

    input:
    path input_mdata

    output:
    path "inference_mudata.h5mu", emit: inference_mudata

    script:
        """
        eval "\$(micromamba shell hook --shell bash)"
        micromamba activate 
        micromamba run -n perturbo_env python /home/jovyan/CRISPR-jamboree_0926/IGVF_Workflows/bin/perturbo_inference.py ${input_mdata} inference_mudata.h5mu
        micromamba deactivate
        """
}
