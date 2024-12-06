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

### Pipeline Setup

```bash
# Clone the repository from GitHub 
git clone https://github.com/jiangsizhu1201/CRISPR_Pipeline.git

# Navigate to pipeline directory
cd IGVF_Workflows
```

### Perturbo Installation

1. Set up environment:
    ```bash
    # Create and activate a new conda environment for Perturbo
    conda create -n perturbo_env
    conda activate perturbo_env

    # Install Perturbo package from wheel file located in `install_perturbo` directory
    pip install perturbo-0.0.1-py3-none-any.whl

    # Exit environment after installation
    conda deactivate
    ```

2. Verify installation:
    ```bash
    # Activate environment and check Perturbo version
    conda activate perturbo_env
    python -c "import perturbo; print(perturbo.__version__)"
    ```

3. Configure workflow:
    ```bash
    # For example:
    ## Update pipeline configuration file with Perturbo environment path
    withName:inference_perturbo {
        conda = '/path/to/perturbo_env'  # Set to your conda environment location
    }
    ```

## Input Requirements

### Required Data Structure

**Essential Directories:**

1. **fastq_files/**: Raw sequencing data
2. **yaml_files/**: Sequence specification files
3. **Metadata files**: Analysis information

This pipeline requires a specific data structure to function properly. Below is an overview of the required directory organization:

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

For detailed specifications, see our [documentation](https://docs.google.com/document/d/1Z1SOlekIE5uGyXW41XxnszxaYdSw0wdAOUVzfy3fj3M/edit?tab=t.0#heading=h.ctbx1w9hj619).

## Running the Pipeline

The configuration looks mostly good, but let's clarify the resource specifications since they're not just for kallisto mapping:

### Configuration

1. Edit `configs/pipeline.config` to specify:
   - Input data paths
   - Analysis parameters

2. Edit `input.config` to specify:
   - Computing resources
   - Container(singularity/docker)/environment(conda)

    ```bash
    # Resource configuration example:
    withName:process_name {
        cpus = 4               # Number of CPU cores per mapping process (default: 4)
        memory = 64.GB         # RAM allocation per mapping process (default: 64GB)
        container = ''         # Container path (Singularity/Docker)
        conda = ''            # Conda environment path
    }
    ```

3. Hardware requirements:
   - GPU: Required for Perturbo inference
   - Adjust resources based on data size

### Execution

```bash
# Set execute permissions
chmod +x bin/*

# Run pipeline with conda environment
nextflow run main.nf -c input.config -with-conda
```

## Output Description

The output files will be generated in the `pipeline_outputs` and `pipeline_dashboard` directory.

### Generated Files

Within the `pipeline_outputs` directory, you will find:

- inference_mudata.h5mu - MuData format output
- per_element_output.tsv - Per-element analysis
- per_guide_output.tsv - Per-guide analysis

**Structure:**

```
pipeline_outputs/
├── inference_mudata.h5mu    
├── per_element_output.tsv 
├── per_guide_output.tsv 
```

For details, see our [documentation](https://docs.google.com/document/d/1Z1SOlekIE5uGyXW41XxnszxaYdSw0wdAOUVzfy3fj3M/edit?tab=t.0#heading=h.ctbx1w9hj619).

### Generated Figures

The pipeline produces several figures:

Within the `pipeline_dashboard` directory, you will find:

1. **Evaluation Output**:
   - `network_plot.png`: Gene interaction networks visualization.
   - `volcano_plot.png`: gRNA-gene pairs analysis.
   - IGV files (`.bedgraph` and `bedpe`): Genome browser visualization files.

2. **Analysis Figures**:
   - `knee_plot_scRNA.png`: Knee plot of UMI counts vs. barcode index.
   - `scatterplot_scrna.png`: Scatterplot of total counts vs. genes detected, colored by mitochondrial content.
   - `violin_plot.png`: Distribution of gene counts, total counts, and mitochondrial content.
   - `scRNA_barcodes_UMI_thresholds.png`: Number of scRNA barcodes using different Total UMI thresholds.
   - `guides_per_cell_histogram.png`: Histogram of guides per cell.
   - `cells_per_guide_histogram.png`: Histogram of cells per guide.
   - `guides_UMI_thresholds.png`: Simulating the final number of cells with assigned guides using different minimal number thresholds (at least one guide > threshold value). (Use it to inspect how many cells would have assigned guides. This can be used to check if the final number of cells with guides fit with your expected number of cells)
   - `guides_UMI_thresholds.png`: Histogram of the number of sgRNA represented per cell
   - `cells_per_htp_barplot.png`: Number of Cells across Different HTOs
   - `umap_hto.png`: UMAP Clustering of Cells Based on HTOs (The dimensions represent the distribution of HTOs in each cell)
   - `umap_hto_singlets.png`: UMAP Clustering of Cells Based on HTOs (multiplets removed)

3. **seqSpec Plots**:

   - `seqSpec_check_plots.png`: The frequency of each nucleotides along the Read 1 (Use to inspect the expected read parts with their expected signature) and Read 2 (Use to inspect the expected read parts with their expected signature).

**Structure:**
```
pipeline_dashboard/
├── dashboard.html                          # Interactive dashboard
├── evaluation_output/                      
│   ├── network_plot.png                   
│   ├── volcano_plot.png                   
│   ├── igv.bedgraph                      
│   └── igv.bedpe                         
│
├── figures/                               
│   ├── knee_plot_scRNA.png               
│   ├── scatterplot_scrna.png             
│   ├── violin_plot.png                   
│   ├── scRNA_barcodes_UMI_thresholds.png 
│   ├── guides_per_cell_histogram.png     
│   ├── cells_per_guide_histogram.png     
│   ├── guides_UMI_thresholds.png         
│   ├── cells_per_htp_barplot.png         
│   ├── umap_hto.png                      
│   └── umap_hto_singlets.png             
│
├── guide_seqSpec_plots/                   
│   └── seqSpec_check_plots.png           
│
└── hashing_seqSpec_plots/                 
    └── seqSpec_check_plots.png           
```


