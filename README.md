# CRISPR-Jamboree

A comprehensive pipeline for single-cell Perturb-Seq analysis, enabling robust testing and demonstration of CRISPR screening data processing at single-cell resolution.

## Getting Started

### Prerequisites

Before running the pipeline, ensure you have the following dependencies installed:

1. **Nextflow** (Workflow Manager)
   ```bash
   conda install bioconda::nextflow
   ```

2. **Mamba** (Package Manager)
   ```bash
   conda install conda-forge::mamba
   ```

3. **Singularity** (Container Platform)
   - Must be available on your execution environment

## Installation Guide

### Repository Setup

```bash
# Clone the repository
git clone https://github.com/jiangsizhu1201/CRISPR_Pipeline.git
cd IGVF_Workflows
```

### Perturbo Installation

1. Create and activate Perturbo environment
```bash
conda create -n perturbo_env
conda activate perturbo_env

# Install from wheel file
pip install perturbo-0.0.1-py3-none-any.whl
conda deactivate
```
```bash
# Verify Perturbo installation
conda activate perturbo_env
python -c "import perturbo; print(perturbo.__version__)"
```

2. Configure environment path in `input.config`:
```bash
withName:inference_perturbo {
    conda = '/path/to/perturbo_env'  # Replace with your actual path
}
```



Would you like me to add troubleshooting steps or common installation issues?

## Input Requirements

### Required Data Structure
```
example_data/
├── fastq_files/                        # Raw sequencing data
│   ├── {sample}_R1.fastq.gz      # Read 1: Cell barcode and UMI
│   └── {sample}_R2.fastq.gz      # Read 2: Transcript sequence
│
├── yaml_files/                      # SeqSpec yaml files
│   ├── seqspec1.yaml             # Read 1 structure
│   └── seqspec2.yaml             # Read 2 structure
│
├── guide_metadata.tsv              #  TSV file contaiing gRNA metadata
├── hash_metadata.tsv             # TSV file contaiing cell hashing metadata (if applicable)
└── pairs_to_test.csv            # CSV file defining perturbation comparisons (if applicable)
```

For detailed file format specifications and examples, please refer to our [documentation](https://docs.google.com/document/d/1Z1SOlekIE5uGyXW41XxnszxaYdSw0wdAOUVzfy3fj3M/edit?tab=t.0#heading=h.ctbx1w9hj619).

## Running the Pipeline

### Configuration

1. Edit `configs/pipeline.config` to specify:
   - Input data paths
   - Analysis parameters

2. Edit `input.config` to specify:
   - computing resources
   - container(singularity/docker)/envrionment(conda)

2. Verify resource settings:
   - Default: 4 CPUs, 64GB RAM
   - GPU required for Perturbo
   - Adjust based on your data size

### Execution

```bash
chmod +x bin/*
nextflow run main.nf -c input.config -with-conda
```

## Output Description

The output files will be generated in the `pipeline_outputs` directory.

Within the `pipeline_outputs` directory, you will find:

- inference_mudata.h5mu - MuData format output
- per_element_output.tsv - Per-element analysis
- per_guide_output.tsv - Per-guide analysis

```
pipeline_outputs/
├── inference_mudata.h5mu    # MuData format output
├── per_element_output.tsv # Per-element analysis
├── per_guide_output.tsv # Per-guide analysis
└── figures/        # Generated visualizations
```

## Generated Figures

The pipeline produces several figures:

1. 

2. 
