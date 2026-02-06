#!/usr/bin/env python3

import os
import pandas as pd
import glob
import argparse
from pathlib import Path

def merge_esviritu_results(input_dir, output_file):
    """
    Merge EsViritu assembly summary files from all samples into a pivot table format
    
    Args:
        input_dir: Directory containing sample subdirectories with EsViritu results
        output_file: Path to output merged TSV file
    """
    
    # Find all assembly summary files
    assembly_files = glob.glob(os.path.join(input_dir, "*", "*.detected_virus.assembly_summary.tsv"))
    
    if not assembly_files:
        print(f"No assembly summary files found in {input_dir}")
        return
    
    print(f"Found {len(assembly_files)} assembly summary files")
    
    # List to store all dataframes
    all_dfs = []
    
    for file_path in assembly_files:
        try:
            # Read the TSV file
            df = pd.read_csv(file_path, sep='\t')
            
            # Extract sample name from file path
            sample_name = os.path.basename(os.path.dirname(file_path))
            
            print(f"Processing {sample_name}: {len(df)} viral detections")
            
            # Add sample name for merging
            df['sample_name'] = sample_name
            
            # Add the dataframe to our list
            all_dfs.append(df)
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue
    
    if not all_dfs:
        print("No valid assembly summary files could be processed")
        return
    
    # Concatenate all dataframes
    all_data = pd.concat(all_dfs, ignore_index=True)
    
    print(f"Total data shape: {all_data.shape}")
    print(f"Total viral detections across all samples: {len(all_data)}")
    print(f"Unique assemblies: {len(all_data['Assembly'].unique())}")
    
    # Create pivot table: Assembly as rows, samples as columns for RPKMF values
    # First, get the metadata columns (taxonomic info, etc.) for each assembly
    metadata_cols = ['Asm_length', 'kingdom', 'phylum', 'tclass', 'order', 
                     'family', 'genus', 'species', 'subspecies', 'Accession', 'Segment']
    
    # Get unique assembly metadata (take first occurrence of each assembly)
    assembly_metadata = all_data.groupby('Assembly')[metadata_cols].first().reset_index()
    
    # Create pivot table for RPKMF values
    rpkmf_pivot = all_data.pivot_table(
        index='Assembly', 
        columns='sample_name', 
        values='RPKMF', 
        fill_value=0
    ).reset_index()
    
    # Merge metadata with RPKMF pivot table
    merged_table = pd.merge(assembly_metadata, rpkmf_pivot, on='Assembly', how='left')
    
    # Reorder columns: metadata first, then sample RPKMF columns
    all_metadata_cols = ['Assembly'] + metadata_cols
    sample_columns = [col for col in merged_table.columns if col not in all_metadata_cols]
    
    # Sort sample columns numerically by extracting the number at the end
    def extract_number(sample_name):
        import re
        # Extract number from end of sample name (e.g., VIR_GYG_1 -> 1, VIR_C_JY_46 -> 46)
        match = re.search(r'_(\d+)$', sample_name)
        if match:
            return int(match.group(1))
        else:
            return 999999  # Put samples without numbers at the end
    
    sample_columns_sorted = sorted(sample_columns, key=extract_number)
    final_columns = all_metadata_cols + sample_columns_sorted
    merged_table = merged_table[final_columns]
    
    print(f"Final merged table shape: {merged_table.shape}")
    print(f"Columns: {list(merged_table.columns)}")
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Save merged table to TSV
    merged_table.to_csv(output_file, sep='\t', index=False)
    
    print(f"Merged results saved to: {output_file}")
    print(f"Table contains {len(merged_table)} unique viral assemblies across {len(sample_columns)} samples")

def main():
    parser = argparse.ArgumentParser(description='Merge EsViritu assembly summary results from all samples')
    parser.add_argument('--input-dir', required=True, help='Directory containing sample subdirectories with EsViritu results')
    parser.add_argument('--output', required=True, help='Output Excel file path')
    
    args = parser.parse_args()
    
    merge_esviritu_results(args.input_dir, args.output)

if __name__ == "__main__":
    main()