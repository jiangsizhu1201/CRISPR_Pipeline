nextflow.enable.dsl=2

include { evaluation_plot } from './processes/evaluation_plot.nf'
include { evaluation_undefined_plot } from './processes/evaluation_undefined_plot.nf'

workflow evaluation_pipeline {

    take:
    gencode_gtf
    inference_mudata

    main:
    if (params.user_central_nodes == 'undefined') {
        Evaluation = evaluation_undefined_plot(inference_mudata, gencode_gtf, params.central_nodes_num)
    }
    else {
        Evaluation = evaluation_plot(inference_mudata, params.user_central_nodes, gencode_gtf)
    }

    emit:
    evaluation_output_dir = Evaluation.evaluation_output
}
