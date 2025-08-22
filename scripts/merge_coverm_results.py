#!/usr/bin/env python3

import os
import pandas as pd
import glob
import argparse
from pathlib import Path

def merge_coverm_results(input_dir, output_dir, metrics=None):
    """
    Merge count, coverage, and TPM results from CoverM analysis
    
    Args:
        input_dir: Directory containing CoverM results (should have count/, coverage/, tpm/ subdirs)
        output_dir: Directory to save merged results
        metrics: List of metrics to process (default: ['count', 'coverage', 'tpm'])
    """
    if metrics is None:
        metrics = ['count', 'coverage', 'tpm']
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize empty dataframes for each metric
    count_df = None
    coverage_df = None
    tpm_df = None
    
    # Get all sample directories from the first available metric
    sample_names = []
    for metric in metrics:
        metric_dirs = glob.glob(os.path.join(input_dir, metric, "*"))
        if metric_dirs:
            sample_names = [os.path.basename(d) for d in metric_dirs]
            break
    
    if not sample_names:
        print(f"Error: No sample directories found in {input_dir}")
        return None
    
    print(f"Found {len(sample_names)} samples to process")
    
    # Process each metric type
    for metric in metrics:
        print(f"Processing {metric} data...")
        combined_df = None
        
        for sample in sample_names:
            # Construct file path
            file_path = os.path.join(input_dir, metric, sample, f"{sample}.{metric}.tsv")
            col_name = f"{sample}_{metric}"
            
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

def main():
    parser = argparse.ArgumentParser(
        description="Merge CoverM count, coverage, and TPM results from multiple samples"
    )
    
    parser.add_argument(
        "-i", "--input-dir",
        required=True,
        help="Input directory containing CoverM results (should have count/, coverage/, tpm/ subdirs)"
    )
    
    parser.add_argument(
        "-o", "--output-dir", 
        required=True,
        help="Output directory to save merged results"
    )
    
    parser.add_argument(
        "-m", "--metrics",
        nargs="+",
        default=["count", "coverage", "tpm"],
        choices=["count", "coverage", "tpm"],
        help="Metrics to process (default: count coverage tpm)"
    )
    
    parser.add_argument(
        "--comprehensive",
        action="store_true",
        default=True,
        help="Create comprehensive table with all metrics (default: True)"
    )
    
    parser.add_argument(
        "--no-comprehensive",
        action="store_true",
        help="Skip creating comprehensive table"
    )
    
    args = parser.parse_args()
    
    # Handle comprehensive table flag
    create_comprehensive = args.comprehensive and not args.no_comprehensive
    
    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Processing metrics: {args.metrics}")
    print(f"Create comprehensive table: {create_comprehensive}")
    print()
    
    result = merge_coverm_results(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        metrics=args.metrics
    )
    
    if result is not None:
        print("\n✅ Merge completed successfully!")
    else:
        print("\n❌ Merge failed!")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())