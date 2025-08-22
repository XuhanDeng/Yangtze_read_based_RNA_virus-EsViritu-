#!/usr/bin/env python3

import os
import pandas as pd
import glob
import argparse
from pathlib import Path

def merge_coverm_results(input_dir, output_dir, metrics=None, reference_csv=None):
    """
    Merge count, coverage, and TPM results from CoverM analysis
    
    Args:
        input_dir: Directory containing CoverM results (should have count/, coverage/, tpm/ subdirs)
        output_dir: Directory to save merged results
        metrics: List of metrics to process (default: ['count', 'coverage', 'tpm'])
        reference_csv: Path to NCBI virus reference CSV file to merge with results
    """
    if metrics is None:
        metrics = ['count', 'coverage', 'tpm']
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Load reference CSV if provided
    reference_df = None
    if reference_csv and os.path.exists(reference_csv):
        try:
            reference_df = pd.read_csv(reference_csv)
            print(f"Loaded reference CSV: {reference_csv}")
            print(f"Reference CSV shape: {reference_df.shape}")
        except Exception as e:
            print(f"Warning: Could not load reference CSV {reference_csv}: {e}")
            reference_df = None
    elif reference_csv:
        print(f"Warning: Reference CSV file not found: {reference_csv}")
    
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
    
    # Sort samples by numeric suffix (extract number at the end of sample name)
    def extract_numeric_suffix(sample_name):
        import re
        # Find the last number in the sample name
        match = re.search(r'(\d+)$', sample_name)
        return int(match.group(1)) if match else float('inf')
    
    sample_names = sorted(sample_names, key=extract_numeric_suffix)
    print(f"First 10 samples in numeric order: {sample_names[:10]}")
    
    print(f"Found {len(sample_names)} samples to process")
    print(f"Sample order: {sample_names[:5]}..." if len(sample_names) > 5 else f"Sample order: {sample_names}")
    
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
        
        # Fill NaN values with 0 and store the combined dataframe
        if combined_df is not None:
            combined_df = combined_df.fillna(0)
        
        if metric == 'count':
            count_df = combined_df
        elif metric == 'coverage':
            coverage_df = combined_df
        else:
            tpm_df = combined_df
    
    # Reorder columns to match sample order for each metric dataframe
    def reorder_columns(df, sample_names, metric):
        if df is None:
            return None
        # Create ordered column list: Contig first, then samples in numeric order
        ordered_cols = ['Contig'] + [f"{sample}_{metric}" for sample in sample_names]
        # Only keep columns that actually exist in the dataframe
        existing_cols = [col for col in ordered_cols if col in df.columns]
        return df[existing_cols]
    
    # Reorder columns in each dataframe
    count_df = reorder_columns(count_df, sample_names, 'count')
    coverage_df = reorder_columns(coverage_df, sample_names, 'coverage') 
    tpm_df = reorder_columns(tpm_df, sample_names, 'tpm')
    
    # Function to merge with reference data
    def merge_with_reference(df, reference_df, suffix=""):
        if df is None or reference_df is None:
            return df
        try:
            # Create a copy of df with version-stripped Contig for matching
            df_copy = df.copy()
            # Strip version numbers from Contig IDs (e.g., NC_000898.1 -> NC_000898)
            df_copy['Contig_base'] = df_copy['Contig'].str.replace(r'\.\d+$', '', regex=True)
            
            # Merge on Contig_base = Accession 
            merged = pd.merge(reference_df, df_copy, left_on='Accession', right_on='Contig_base', how='right')
            
            # Keep original Contig column and drop the temporary Contig_base
            if 'Contig_base' in merged.columns:
                merged = merged.drop('Contig_base', axis=1)
            
            # Count successful matches
            matches = merged['Accession'].notna().sum()
            total = len(merged)
            print(f"  Merged with reference data{suffix}: {merged.shape} ({matches}/{total} matches)")
            return merged
        except Exception as e:
            print(f"  Warning: Could not merge with reference data{suffix}: {e}")
            return df

    # Save individual metric files (with reference data if available)
    for df, metric in [(count_df, 'count'), (coverage_df, 'coverage'), (tpm_df, 'tpm')]:
        if df is not None:
            # Merge with reference data if available
            if reference_df is not None:
                df_with_ref = merge_with_reference(df, reference_df, f" for {metric}")
                ref_output_file = os.path.join(output_dir, f"merged_{metric}_with_reference_table.tsv")
                df_with_ref.to_csv(ref_output_file, sep='\t', index=False)
                print(f"Saved {metric} table with reference: {ref_output_file}")
                print(f"  Shape: {df_with_ref.shape}")
            
            # Also save original without reference
            output_file = os.path.join(output_dir, f"merged_{metric}_table.tsv")
            df.to_csv(output_file, sep='\t', index=False)
            print(f"Saved {metric} table: {output_file}")
            print(f"  Shape: {df.shape}")
    
    # Create a comprehensive merged table with all metrics
    if count_df is not None and coverage_df is not None and tpm_df is not None:
        print("Creating comprehensive merged table...")
        
        # Start with count data
        merged_all = count_df.copy()
        
        # Add coverage data (maintaining sample order)
        coverage_cols = [col for col in coverage_df.columns if col != 'Contig']
        for col in coverage_cols:
            merged_all = pd.merge(merged_all, coverage_df[['Contig', col]], on='Contig', how='outer')
        
        # Add TPM data (maintaining sample order)
        tpm_cols = [col for col in tpm_df.columns if col != 'Contig']
        for col in tpm_cols:
            merged_all = pd.merge(merged_all, tpm_df[['Contig', col]], on='Contig', how='outer')
        
        # Reorder columns in comprehensive table: Contig, then all count, all coverage, all tpm
        ordered_comprehensive_cols = ['Contig']
        for sample in sample_names:
            # Add count, coverage, tpm for each sample in order
            for metric in ['count', 'coverage', 'tpm']:
                col_name = f"{sample}_{metric}"
                if col_name in merged_all.columns:
                    ordered_comprehensive_cols.append(col_name)
        
        # Reorder the comprehensive dataframe
        merged_all = merged_all[ordered_comprehensive_cols]
        
        # Fill NaN values with 0
        merged_all = merged_all.fillna(0)
        
        # Save comprehensive table with reference data if available
        if reference_df is not None:
            comprehensive_with_ref = merge_with_reference(merged_all, reference_df, " for comprehensive")
            comprehensive_ref_output = os.path.join(output_dir, "merged_comprehensive_with_reference_table.tsv")
            comprehensive_with_ref.to_csv(comprehensive_ref_output, sep='\t', index=False)
            print(f"Saved comprehensive table with reference: {comprehensive_ref_output}")
            print(f"  Shape: {comprehensive_with_ref.shape}")
        
        # Save original comprehensive table without reference
        comprehensive_output = os.path.join(output_dir, "merged_comprehensive_table.tsv")
        merged_all.to_csv(comprehensive_output, sep='\t', index=False)
        print(f"Saved comprehensive table: {comprehensive_output}")
        print(f"  Shape: {merged_all.shape}")
        
        # Display summary statistics
        print("\nSummary:")
        print(f"Total contigs: {len(merged_all)}")
        print(f"Total samples: {len(sample_names)}")
        print(f"Columns in comprehensive table: {len(merged_all.columns)}")
        if reference_df is not None:
            print(f"Reference data columns: {len(reference_df.columns)}")
        
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
    
    parser.add_argument(
        "-r", "--reference-csv",
        help="Path to NCBI virus reference CSV file to merge with results"
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
        metrics=args.metrics,
        reference_csv=args.reference_csv
    )
    
    if result is not None:
        print("\n✅ Merge completed successfully!")
    else:
        print("\n❌ Merge failed!")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())