
nextflow.enable.dsl=2

include { PreprocessAnnData } from './processes/PreprocessAnnData.nf'
include { CreateMuData_HASHING } from './processes/CreateMuData_HASHING.nf' 
include { demultiplex } from './processes/demultiplex.nf'
include { filter_hashing } from './processes/filter_hashing.nf'
include { hashing_concat } from './processes/hashing_concat.nf'
include { doublets_scrub_HASHING } from './processes/doublets_scrub_HASHING.nf'
include { guide_assignment_cleanser } from './processes/guide_assignment_cleanser.nf'
include { guide_assignment_sceptre } from './processes/guide_assignment_sceptre.nf'
include { guide_assignment_mudata } from './processes/guide_assignment_mudata.nf'
include { downloadGTF } from './processes/downloadGTF.nf'
include { prepare_guide_inference } from './processes/prepare_guide_inference.nf'
include { prepare_user_guide_inference } from './processes/prepare_user_guide_inference.nf'
include { inference_sceptre } from './processes/inference_sceptre.nf'
include { inference_perturbo } from './processes/inference_perturbo.nf'
include { inference_mudata } from './processes/inference_mudata.nf'
include { mergedResults } from './processes/mergedResults.nf'

workflow process_mudata_pipeline_HASHING {

    take:
    concat_anndata_rna
    trans_out_dir
    concat_anndata_guide
    guide_out_dir
    concat_anndata_hashing
    hashing_out_dir
    covariate_string

    main:

    Preprocessed_AnnData = PreprocessAnnData(
        concat_anndata_rna,
        trans_out_dir.flatten().first(),
        params.min_genes,
        params.min_cells,
        params.pct_mito,
        params.transcriptome
        )

    Hashing_Filtered = filter_hashing(
        Preprocessed_AnnData.filtered_anndata_rna, 
        concat_anndata_hashing
        )

    hashing_filtered_anndata_collected = Hashing_Filtered.hashing_filtered_anndata.collect()
    hashing_filtered_anndata_collected.view()

    Demultiplex = demultiplex(hashing_filtered_anndata_collected)

    hashing_demux_anndata_collected =Demultiplex.hashing_demux_anndata.collect()
    hashing_demux_anndata_collected.view()

    Hashing_Concat = hashing_concat(hashing_demux_anndata_collected)
    
    GTF_Reference = downloadGTF(params.gtf_url)

    MuData = CreateMuData_HASHING(
        Preprocessed_AnnData.filtered_anndata_rna,
        concat_anndata_guide, 
        Hashing_Concat.concatenated_hashing_demux,
        file(params.guide_metadata),
        GTF_Reference.gencode_gtf,
        params.moi
        )

    MuData_Doublets = doublets_scrub_HASHING(MuData.mudata) 
    
    if (params.assignment_method == "cleanser") {
    Guide_Assignment = guide_assignment_cleanser(MuData_Doublets.mudata_doublet, params.THRESHOLD)}
    else if (params.assignment_method == "sceptre") {
    Guide_Assignment_Matrix = guide_assignment_sceptre(MuData_Doublets.mudata_doublet)
    Guide_Assignment = guide_assignment_mudata(
        Guide_Assignment_Matrix.guide_assignment_matrix, 
        MuData_Doublets.mudata_doublet)
    }

    if (params.inference_option == 'predefined_pairs') {
        PrepareInference = prepare_user_guide_inference(
            Guide_Assignment.guide_assignment_mudata_output,
            file(params.user_inference)
        )}
    else {
        PrepareInference = prepare_guide_inference(
            Guide_Assignment.guide_assignment_mudata_output,
            GTF_Reference.gencode_gtf,
            params.distance_from_center
        )}

    if (params.inference_method == "sceptre"){
        TestResults = inference_sceptre(PrepareInference.mudata_inference_input, covariate_string)
        GuideInference = inference_mudata(TestResults.test_results, PrepareInference.mudata_inference_input, params.inference_method)
    }
    else if (params.inference_method == "perturbo"){
        GuideInference = inference_perturbo(PrepareInference.mudata_inference_input, params.inference_method)
    }
    else if (params.inference_method == "sceptre,perturbo") {
        SceptreResults = inference_sceptre(PrepareInference.mudata_inference_input, covariate_string)
        PerturboResults = inference_perturbo(PrepareInference.mudata_inference_input,  "perturbo")
        GuideInference = mergedResults(SceptreResults.test_results, PerturboResults.inference_mudata)
    }

    emit:
    inference_mudata = GuideInference.inference_mudata
    gencode_gtf = GTF_Reference.gencode_gtf
    figures_dir = Preprocessed_AnnData.figures_dir
    adata_rna = Preprocessed_AnnData.adata_rna
    filtered_anndata_rna = Preprocessed_AnnData.filtered_anndata_rna
    adata_guide = MuData.adata_guide
    adata_hashing = Hashing_Filtered.adata_hashing
    adata_demux = Hashing_Concat.concatenated_hashing_demux
}
