nextflow.enable.dsl=2

include { skipGenomeDownload } from './processes/skipGenomeDownload.nf'
include { downloadGenome } from './processes/downloadGenome.nf'
include { prepare_covariate } from './processes/prepare_covariate.nf'

workflow prepare_mapping_pipeline {
    if (file(params.genome_local_path).exists()) {
        Genome = skipGenomeDownload(file(params.genome_local_path))
    }
    else {
        Genome = downloadGenome(params.genome_download_path)
    }

    Prepare_covariate = prepare_covariate(params.covariate_list)

    emit:
    genome = Genome.genome
    parsed_covariate_file =  Prepare_covariate.parsed_covariate_file
    covariate_string =  Prepare_covariate.covariate_string
}
