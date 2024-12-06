nextflow.enable.dsl=2

include { seqSpecParser } from './processes/seqSpecParser.nf'
include { downloadReference } from './processes/downloadReference.nf'
include { mappingscRNA } from './processes/mappingscRNA.nf'
include { anndata_concat } from './processes/anndata_concat.nf'

workflow mapping_rna_pipeline{
    take:
    parsed_covariate_file

    main:
    SeqSpecResult = seqSpecParser(
        file("${params.seqspecs_directory}/${params.scRNA_seqspec_yaml}"),
        file(params.seqspecs_directory),
        'rna'
        )

    DownloadRefResult = downloadReference(params.transcriptome)

    fastq_files_ch = Channel.fromPath(params.fastq_files_rna)
    batch_ch = Channel.fromList(params.batch)

    batch_fastq_ch = batch_ch.merge(fastq_files_ch)

    MappingOut = mappingscRNA(
        batch_fastq_ch,
        DownloadRefResult.transcriptome_idx,
        DownloadRefResult.t2g_transcriptome_index,
        SeqSpecResult.parsed_seqspec,
        SeqSpecResult.barcode_file
        )
    
    ks_transcripts_out_dir_collected = MappingOut.ks_transcripts_out_dir.collect()
    ks_transcripts_out_dir_collected.view()

    AnndataConcatenate = anndata_concat(
        parsed_covariate_file,
        ks_transcripts_out_dir_collected
    )

    emit:
    trans_out_dir = MappingOut.ks_transcripts_out_dir
    ks_transcripts_out_dir_collected = MappingOut.ks_transcripts_out_dir.collect()
    concat_anndata_rna = AnndataConcatenate.concat_anndata
}
