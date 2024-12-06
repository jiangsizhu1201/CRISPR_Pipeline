nextflow.enable.dsl=2

include { seqSpecCheck_pipeline } from './seqSpecCheck_pipeline.nf'
include { seqSpecCheck_pipeline_HASHING } from './seqSpecCheck_pipeline_HASHING.nf'
include { prepare_mapping_pipeline } from './prepare_mapping_pipeline.nf'
include { mapping_rna_pipeline } from './mapping_rna_pipeline.nf'
include { mapping_guide_pipeline } from './mapping_guide_pipeline.nf'
include { mapping_hashing_pipeline } from './mapping_hashing_pipeline.nf'
include { process_mudata_pipeline_HASHING } from './process_mudata_pipeline_HASHING.nf'
include { process_mudata_pipeline } from './process_mudata_pipeline.nf'
include { evaluation_pipeline } from './evaluation_pipeline.nf'
include { dashboard_pipeline_HASHING } from './dashboard_pipeline_HASHING.nf'
include { dashboard_pipeline } from './dashboard_pipeline.nf'

workflow {
  
  if (params.DATASET_HASHING == "true"){
    seqSpecCheck_pipeline_HASHING() 
    }
  else {
    seqSpecCheck_pipeline()
    }

  prepare_mapping_pipeline()

  mapping_rna_pipeline(
    prepare_mapping_pipeline.out.parsed_covariate_file
    )
  mapping_guide_pipeline(
    prepare_mapping_pipeline.out.parsed_covariate_file,
    prepare_mapping_pipeline.out.genome
    )

  if (params.DATASET_HASHING == "true"){

    mapping_hashing_pipeline(
      prepare_mapping_pipeline.out.parsed_covariate_file,
      prepare_mapping_pipeline.out.genome
      )

    process_mudata_pipeline_HASHING(
      mapping_rna_pipeline.out.concat_anndata_rna,
      mapping_rna_pipeline.out.trans_out_dir,
      mapping_guide_pipeline.out.concat_anndata_guide,
      mapping_guide_pipeline.out.guide_out_dir,
      mapping_hashing_pipeline.out.concat_anndata_hashing,
      mapping_hashing_pipeline.out.hashing_out_dir,
      prepare_mapping_pipeline.out.covariate_string
    )

    evaluation_pipeline (
      process_mudata_pipeline_HASHING.out.gencode_gtf,
      process_mudata_pipeline_HASHING.out.inference_mudata
    )

    dashboard_pipeline_HASHING (
      seqSpecCheck_pipeline_HASHING.out.guide_seqSpecCheck_plots,
      seqSpecCheck_pipeline_HASHING.out.guide_position_table,
      seqSpecCheck_pipeline_HASHING.out.hashing_seqSpecCheck_plots,
      seqSpecCheck_pipeline_HASHING.out.hashing_position_table,
      process_mudata_pipeline_HASHING.out.adata_rna,
      process_mudata_pipeline_HASHING.out.filtered_anndata_rna,
      mapping_rna_pipeline.out.ks_transcripts_out_dir_collected,
      process_mudata_pipeline_HASHING.out.adata_guide,
      mapping_guide_pipeline.out.ks_guide_out_dir_collected,
      process_mudata_pipeline_HASHING.out.adata_hashing,
      mapping_hashing_pipeline.out.ks_hashing_out_dir_collected,
      process_mudata_pipeline_HASHING.out.adata_demux,
      process_mudata_pipeline_HASHING.out.inference_mudata,
      process_mudata_pipeline_HASHING.out.figures_dir,
      evaluation_pipeline.out.evaluation_output_dir
    )
    
  }
  else {
    process_mudata_pipeline(
      mapping_rna_pipeline.out.concat_anndata_rna,
      mapping_rna_pipeline.out.trans_out_dir,
      mapping_guide_pipeline.out.concat_anndata_guide,
      mapping_guide_pipeline.out.guide_out_dir,
      prepare_mapping_pipeline.out.covariate_string
    )

    evaluation_pipeline (
      process_mudata_pipeline.out.gencode_gtf,
      process_mudata_pipeline.out.inference_mudata
    )

    dashboard_pipeline (
      seqSpecCheck_pipeline.out.guide_seqSpecCheck_plots,
      seqSpecCheck_pipeline.out.guide_position_table,
      process_mudata_pipeline.out.adata_rna,
      process_mudata_pipeline.out.filtered_anndata_rna,
      mapping_rna_pipeline.out.ks_transcripts_out_dir_collected,
      process_mudata_pipeline.out.adata_guide,
      mapping_guide_pipeline.out.ks_guide_out_dir_collected,
      process_mudata_pipeline.out.inference_mudata,
      process_mudata_pipeline.out.figures_dir,
      evaluation_pipeline.out.evaluation_output_dir
    )
  }

}
