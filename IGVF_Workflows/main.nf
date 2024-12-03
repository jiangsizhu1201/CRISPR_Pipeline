nextflow.enable.dsl=2

include { seqSpecCheck_pipeline } from './seqSpecCheck_pipeline.nf'
include { seqSpecCheck_pipeline_HASHING } from './seqSpecCheck_pipeline_HASHING.nf'
include { prepare_mapping_pipeline } from './prepare_mapping_pipeline.nf'
include { mapping_rna_pipeline } from './mapping_rna_pipeline.nf'
include { mapping_guide_pipeline } from './mapping_guide_pipeline.nf'
include { mapping_hashing_pipeline } from './mapping_hashing_pipeline.nf'
include { process_mudata_pipeline_HASHING } from './process_mudata_pipeline_HASHING.nf'
include { process_mudata_pipeline } from './process_mudata_pipeline.nf'

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
    
  }
  else {
    process_mudata_pipeline(
      mapping_rna_pipeline.out.concat_anndata_rna,
      mapping_rna_pipeline.out.trans_out_dir,
      mapping_guide_pipeline.out.concat_anndata_guide,
      mapping_guide_pipeline.out.guide_out_dir,
      prepare_mapping_pipeline.out.covariate_string
    )
  }

}
