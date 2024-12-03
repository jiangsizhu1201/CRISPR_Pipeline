#!/usr/bin/env python3

import sys
import pandas as pd

def parse_tsv_and_generate_params(tsv_file):
    data = pd.read_csv(tsv_file, sep='\t')
    
    # Set DATASET_HASHING 
    has_hash = 'hash' in data['file_modality'].values
    
    params = {
        'DATASET_HASHING': "true" if has_hash else "false",
        'seqspecs_directory': "example_data/yaml_files",
        'user_inference': "example_data/user_pairs_to_test.csv",
        'guide_metadata': "example_data/guide_metadata.tsv",
        'hashing_metadata': "example_data/hash_metadata.tsv",
        'transcriptome': "human",
        'genome_download_path': "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        'genome_local_path': "",
        'gtf_url': "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz",
        'scRNA_seqspec_yaml': "",
        'Guides_seqspec_yaml': "",
        'Hash_seqspec_yaml': "",
        'assignment_method': "sceptre",
        'THRESHOLD': "1",
        'inference_method': "sceptre,perturbo",
        'moi': "undecided",
        'side': "both",
        'grna_integration_strategy': "union",
        'resampling_approximation': "skew_normal",
        'control_group': "default",
        'resampling_mechanism': "default",
        'formula_object': "default",
        'inference_option': "predefined_pairs",
        'distance_from_center': "1000000",
        'min_genes': "",
        'min_cells': "",
        'pct_mito': "",
        'user_central_nodes': "undefined",
        'central_nodes_num': "2"
    }

    # Group by measurement_sets and file_modality
    grouped = data.groupby(['measurement_sets', 'file_modality'])

    rna_files = []
    guide_files = []
    hashing_files = []

    # Process each group
    for (measurement_set, modality), group in grouped:
        file_pairs = [
            f'example_data/fastq_files/{r1}.fastq.gz example_data/fastq_files/{r2}.fastq.gz'
            for r1, r2 in zip(group['R1_path'], group['R2_path'])
        ]
        
        if modality == 'scRNA':
            rna_files.extend(file_pairs)
        elif modality == 'gRNA':
            guide_files.extend(file_pairs)
        elif modality == 'hash':
            hashing_files.extend(file_pairs)

    # Get first batch from gRNA and hash files
    first_guide = guide_files[0].split() if guide_files else ['', '']
    first_hash = hashing_files[0].split() if hashing_files else ['', '']

    # Add file lists to params
    params.update({
        'fastq_files_rna': rna_files,
        'fastq_files_guide': guide_files,
        'fastq_files_hashing': hashing_files,
        'test_guide_fastq_r1': [first_guide[0]],
        'test_guide_fastq_r2': [first_guide[1]],
        'test_hashing_fastq_r1': [first_hash[0]],
        'test_hashing_fastq_r2': [first_hash[1]],
    })

    measurement_to_lane = {}
    for measurement_set in data['measurement_sets'].unique():
        lane = data[data['measurement_sets'] == measurement_set]['Lane'].iloc[0]
        measurement_to_lane[measurement_set] = f"lane{lane}"

    batch_names = []
    lane_values = []
    for i, measurement_set in enumerate(data['measurement_sets'].unique()):
        batch_names.append(f"batch_{chr(97 + i)}")  
        lane_values.append(measurement_to_lane[measurement_set])

    params['batch'] = batch_names

    covariate_list = {
        'batch': batch_names,
        'cov1': lane_values
    }

    return params, covariate_list

def format_output(params, covariate_list):
    output = ["params {"]
    
    # Format params
    for key, value in params.items():
        if isinstance(value, list):
            output.append(f"    {key} = [")
            for item in value:
                output.append(f'        "{item}",')
            output.append("    ]")
        else:
            output.append(f'    {key} = "{value}"')
    output.append("}")
    
    # Format covariate list
    output.append("\nparams.covariate_list = [")
    for key, value in covariate_list.items():
        output.append(f"    {key}: {value},")
    output.append("]")
    
    return "\n".join(output)

def save_config(formatted_output, output_file='pipeline_input.config'):
    with open(output_file, 'w') as f:
        f.write(formatted_output)

def main():
    if len(sys.argv) != 2:
        print("Usage: ./generate_config.py <tsv_file>")
        sys.exit(1)
        
    tsv_file = sys.argv[1]
    
    try:
        params, covariate_list = parse_tsv_and_generate_params(tsv_file)
        formatted_output = format_output(params, covariate_list)
        save_config(formatted_output)
        print(f"Configuration has been successfully saved to pipeline_input.config")
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
