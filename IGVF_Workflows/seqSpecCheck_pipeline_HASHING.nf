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

workflow hashingWorkflow {
    hashing_fastq_r1_ch = Channel.fromPath(params.test_hashing_fastq_r1).collect()
    hashing_fastq_r2_ch = Channel.fromPath(params.test_hashing_fastq_r2).collect()

    hashing_seqSpecCheck = seqSpecCheck(
        hashing_fastq_r1_ch,
        hashing_fastq_r2_ch,
        file(params.hashing_metadata),
        'hashing'
    )

    emit:
    hashing_seqSpecCheck_plots = hashing_seqSpecCheck.seqSpecCheck_plots
    hashing_position_table = hashing_seqSpecCheck.position_table
}

workflow seqSpecCheck_pipeline_HASHING {
    main:
    guide = guideWorkflow()
    hashing = hashingWorkflow()

    emit:
    guide_seqSpecCheck_plots = guide.guide_seqSpecCheck_plots
    guide_position_table = guide.guide_position_table
    hashing_seqSpecCheck_plots = hashing.hashing_seqSpecCheck_plots
    hashing_position_table = hashing.hashing_position_table
}
