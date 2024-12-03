
process inference_perturbo {
    cache 'lenient'
    publishDir './pipeline_outputs'

    input:
    path mudata
    val inference_method

    output:
    path "inference_mudata.h5mu", emit: inference_mudata
    path "per_element_output.tsv", emit: per_element_output
    path "per_guide_output.tsv", emit: per_guide_output

    script:
        def nf_bin = "${projectDir}/bin"
        """
        # Use absolute paths to everything
        ${task.conda}/bin/python ${nf_bin}/perturbo_inference.py ${mudata} inference_mudata.h5mu
        ${task.conda}/bin/python ${nf_bin}/export_output_single.py --mudata inference_mudata.h5mu --inference_method ${inference_method}
        """
}
