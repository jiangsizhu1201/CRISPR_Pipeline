process inference_mudata {
  
    cache 'lenient'
    
    input:
    path test_result
    path mudata

    output:
        path "inference_mudata.h5mu", emit: inference_mudata

    script:
        """
          add_guide_inference.py --test_results_csv ${test_result} --mudata ${mudata}
        """

}
