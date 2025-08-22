#!/usr/bin/env python3

import os
import pandas as pd
import glob
from pathlib import Path

def merge_coverm_results():
    """
    Merge count, coverage, and TPM results from CoverM analysis
    """
    base_dir = "/Users/dengxuhan/Desktop/yangtze/RNA/read_based/results/6_bowtie2/reduntant/3_coverm"
    output_dir = "/Users/dengxuhan/Desktop/yangtze/RNA/read_based/results/6_bowtie2/redundant"
    
    # Initialize empty dataframes for each metric
    count_df = None
    coverage_df = None
    tpm_df = None
    
    # Get all sample directories
    count_dirs = glob.glob(os.path.join(base_dir, "count", "*"))
    sample_names = [os.path.basename(d) for d in count_dirs]
    
    print(f"Found {len(sample_names)} samples to process")
    
    # Process each metric type
    for metric in ['count', 'coverage', 'tpm']:
        print(f"Processing {metric} data...")
        combined_df = None
        
        for sample in sample_names:
            # Construct file path
            if metric == 'count':
                file_path = os.path.join(base_dir, metric, sample, f"{sample}.count.tsv")
                col_name = f"{sample}_count"
            elif metric == 'coverage':
                file_path = os.path.join(base_dir, metric, sample, f"{sample}.coverage.tsv")
                col_name = f"{sample}_coverage"
            else:  # tpm
                file_path = os.path.join(base_dir, metric, sample, f"{sample}.tpm.tsv")
                col_name = f"{sample}_tpm"
            
            # Check if file exists
            if not os.path.exists(file_path):
                print(f"Warning: File not found: {file_path}")
                continue
            
            # Read the TSV file
            try:
                df = pd.read_csv(file_path, sep='\t')
                
                # Rename the second column to sample-specific name
                df.columns = ['Contig', col_name]
                
                # Merge with combined dataframe
                if combined_df is None:
                    combined_df = df
                else:
                    combined_df = pd.merge(combined_df, df, on='Contig', how='outer')
                    
            except Exception as e:
                print(f"Error processing {file_path}: {e}")
                continue
        
        # Store the combined dataframe
        if metric == 'count':
            count_df = combined_df
        elif metric == 'coverage':
            coverage_df = combined_df
        else:
            tpm_df = combined_df
    
    # Save individual metric files
    for df, metric in [(count_df, 'count'), (coverage_df, 'coverage'), (tpm_df, 'tpm')]:
        if df is not None:
            output_file = os.path.join(output_dir, f"merged_{metric}_table.tsv")
            df.to_csv(output_file, sep='\t', index=False)
            print(f"Saved {metric} table: {output_file}")
            print(f"  Shape: {df.shape}")
    
    # Create a comprehensive merged table with all metrics
    if count_df is not None and coverage_df is not None and tpm_df is not None:
        print("Creating comprehensive merged table...")
        
        # Start with count data
        merged_all = count_df.copy()
        
        # Add coverage data
        coverage_cols = [col for col in coverage_df.columns if col != 'Contig']
        for col in coverage_cols:
            merged_all = pd.merge(merged_all, coverage_df[['Contig', col]], on='Contig', how='outer')
        
        # Add TPM data
        tpm_cols = [col for col in tpm_df.columns if col != 'Contig']
        for col in tpm_cols:
            merged_all = pd.merge(merged_all, tpm_df[['Contig', col]], on='Contig', how='outer')
        
        # Fill NaN values with 0
        merged_all = merged_all.fillna(0)
        
        # Save comprehensive table
        comprehensive_output = os.path.join(output_dir, "merged_comprehensive_table.tsv")
        merged_all.to_csv(comprehensive_output, sep='\t', index=False)
        print(f"Saved comprehensive table: {comprehensive_output}")
        print(f"  Shape: {merged_all.shape}")
        
        # Display summary statistics
        print("\nSummary:")
        print(f"Total contigs: {len(merged_all)}")
        print(f"Total samples: {len(sample_names)}")
        print(f"Columns in comprehensive table: {len(merged_all.columns)}")
        
        return merged_all
    
    else:
        print("Error: Could not create comprehensive table - some metric dataframes are missing")
        return None

if __name__ == "__main__":
    result = merge_coverm_results()