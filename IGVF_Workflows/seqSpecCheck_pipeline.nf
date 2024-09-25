process seqSpecCheck {
    cache 'lenient'
    debug true
    
    input:
    path(test_fastq_files_r1)
    path(test_fastq_files_r2)
    path(metadata)
    val data_type
    
    output:
    path "*_seqSpec_plots", emit: seqSpecCheck_plots
    path "*_position_table.csv", emit: position_table
    
    script:
    def r1_files = test_fastq_files_r1.join(' ')
    def r2_files = test_fastq_files_r2.join(' ')
    """
    echo "Checking fastq files for ${data_type}"
    seqSpecCheck.py --read1 ${r1_files} --read2 ${r2_files} --max_reads 100000 --metadata ${metadata} --plot
    mv seqSpec_plots ${data_type}_seqSpec_plots
    mv position_table.csv ${data_type}_position_table.csv
    """
}

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

workflow seqSpecCheck_pipeline {
    main:
    guide = guideWorkflow()
    hashing = hashingWorkflow()

    emit:
    guide_seqSpecCheck_plots = guide.guide_seqSpecCheck_plots
    guide_position_table = guide.guide_position_table
    hashing_seqSpecCheck_plots = hashing.hashing_seqSpecCheck_plots
    hashing_position_table = hashing.hashing_position_table
}
