nextflow.enable.dsl=2

include { seqSpecCheck } from './processes/seqSpecCheck.nf'

workflow guideWorkflow {
    guide_fastq_r1_ch = Channel.fromPath(params.test_guide_fastq_r1).collect()
    guide_fastq_r2_ch = Channel.fromPath(params.test_guide_fastq_r2).collect()
    
    guide_seqSpecCheck = seqSpecCheck(
        guide_fastq_r1_ch,
        guide_fastq_r2_ch,
        file(params.guide_metadata),
        'guide'
    )

    emit:
    guide_seqSpecCheck_plots = guide_seqSpecCheck.seqSpecCheck_plots
    guide_position_table = guide_seqSpecCheck.position_table
}


workflow seqSpecCheck_pipeline {
    main:
    guide = guideWorkflow()

    emit:
    guide_seqSpecCheck_plots = guide.guide_seqSpecCheck_plots
    guide_position_table = guide.guide_position_table
}
