# Pipeline Configuration Documentation

> **IMPORTANT**: All parameters defined in this configuration must have valid, non-empty values. The Nextflow pipeline requires every parameter to be properly set - empty or missing parameters will cause pipeline execution failures. Please ensure all parameters are specified before running the pipeline.

This document explains the configuration parameters for the single-cell RNA sequencing and CRISPR screening analysis pipeline.

## Input Data Parameters

### File Paths
- `user_inference`: Path to CSV file containing user-defined cell pairs for testing
- `guide_metadata`: Path to TSV file containing guide RNA metadata
- `hashing_metadata`: Path to TSV file containing cell hashing metadata

### Dataset Configuration
- `DATASET_HASHING`: Boolean flag to enable/disable dataset hashing ('true'/'false')
- `transcriptome`: Species specification for transcriptome analysis (set to 'human')
- `seqspecs_directory`: Directory containing YAML specification files
- `genome_download_path`: URL for downloading the human genome (hg38) reference
- `genome_local_path`: Local path to store the downloaded genome file
- `gtf_url`: URL for downloading GENCODE gene annotations (v46)

### Sequence Specification Files
- `scRNA_seqspec_yaml`: YAML file specifying single-cell RNA sequencing parameters
- `Guides_seqspec_yaml`: YAML file specifying guide RNA sequencing parameters
- `Hash_seqspec_yaml`: YAML file specifying cell hashing parameters

## Analysis Parameters

### Inference Configuration
- `assignment_method`: Method for guide-cell assignment (set to 'sceptre')
- `THRESHOLD`: Numerical threshold for assignments (set to 1)
- `inference_method`: Methods for statistical inference ('sceptre,perturbo')
- `moi`: Multiplicity of infection setting ('undecided')
- `side`: Direction of statistical tests ('both')
- `grna_integration_strategy`: Strategy for combining guide information ('union')
- `resampling_approximation`: Distribution for resampling ('skew_normal')
- `control_group`: Control group specification ('default')
- `resampling_mechanism`: Mechanism for resampling ('default')
- `formula_object`: Statistical formula specification ('default')

### Analysis Options
- `inference_option`: Type of inference analysis ('predefined_pairs')
- `distance_from_center`: Distance threshold for spatial analysis (1,000,000 base pairs)
- `min_genes`: Minimum number of genes required per cell (500)
- `min_cells`: Minimum number of cells required (3)
- `pct_mito`: Maximum percentage of mitochondrial genes allowed (20%)
- `user_central_nodes`: User-defined central nodes ('undefined')
- `central_nodes_num`: Number of central nodes to consider (2)

## Input Files Configuration

### FASTQ Files
The pipeline accepts three types of FASTQ files:

1. RNA Sequencing Files (`fastq_files_rna`):
   - Paired-end reads organized in pairs
   - Example format: `IGVFFI1946LEGM.fastq.gz, IGVFFI5195OGCL.fastq.gz`

2. Guide RNA Sequencing Files (`fastq_files_guide`):
   - Paired-end reads for guide RNA sequencing
   - Example format: `IGVFFI7706SWGW.fastq.gz, IGVFFI7788FDIR.fastq.gz`

3. Cell Hashing Files (`fastq_files_hashing`):
   - Paired-end reads for cell hashing
   - Example format: `IGVFFI5460OSRQ.fastq.gz, IGVFFI1587BLSX.fastq.gz`

### Test Files
Specific test files are configured for validation:
- Guide RNA test files (R1 and R2)
- Cell hashing test files (R1 and R2)

## Batch Information

### Batch Configuration
- Batch identifiers: 'batch_a', 'batch_b'
- Covariate list includes:
  - Batch information
  - Additional covariates (e.g., lane information)

## Usage Notes

1. All file paths should be relative to the project directory (`${projectDir}`)
2. FASTQ files should be properly paired and organized in the specified directory
3. The pipeline supports both local and remote genome references
4. Ensure all YAML specification files are present in the specified directory
5. Adjust threshold values and parameters according to your experimental design
6. Validate that all parameters have proper values - the pipeline will fail if any parameter is empty or missing
7. Use placeholder values like 'default', 'undefined', or 'undecided' when actual values are not yet determined, but never leave parameters empty

## Parameter Validation Requirements

- Every parameter in the configuration must be explicitly defined
- Use appropriate default values when actual values are not yet determined
- Ensure all file paths point to existing files or locations
- Boolean parameters must be explicitly set to 'true' or 'false'