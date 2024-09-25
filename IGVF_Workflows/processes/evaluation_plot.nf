
process evaluation_plot {
  
    cache 'lenient'
    input:

    path mdata
    val user_central_nodes
    path gencode_gtf

    output:
    path "evaluation_output" , emit: evaluation_output

    script:
            """
            network_plot.py ${mdata} --central_node ${user_central_nodes} 
            volcano_plot.py ${mdata} --log2_fc=1 --p_value=0.05
            igv.py ${mdata} --gtf ${gencode_gtf}
            """
}
