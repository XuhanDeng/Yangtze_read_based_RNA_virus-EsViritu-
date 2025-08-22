#!/usr/bin/env python3

import os
import pandas as pd
import argparse
import glob

def convert_count_to_rpkmf(count_table_path, count_with_ref_path, filtered_reads_dir, 
                          contig_lengths_path, output_rpkmf, output_rpkmf_ref):
    """
    Convert CoverM read counts to RPKMF using EsViritu formula:
    RPKMF = (read_count / (contig_length / 1000)) / (filtered_reads / 1e6)
    """
    
    print("Loading input data...")
    
    # 1. Load contig lengths
    contig_lengths = pd.read_csv(contig_lengths_path, sep='\t')
    print(f"Loaded contig lengths: {contig_lengths.shape}")
    
    # 2. Load filtered read counts for each sample
    filtered_reads = {}
    count_files = glob.glob(os.path.join(filtered_reads_dir, "*_filtered_count.txt"))
    
    for count_file in count_files:
        sample_name = os.path.basename(count_file).replace("_filtered_count.txt", "")
        with open(count_file, 'r') as f:
            count = int(f.read().strip())
            filtered_reads[sample_name] = count
    
    print(f"Loaded filtered read counts for {len(filtered_reads)} samples")
    
    # 3. Process count tables
    def convert_table_to_rpkmf(count_table_path, output_path, table_name):
        print(f"\nProcessing {table_name}...")
        
        # Load count table
        count_df = pd.read_csv(count_table_path, sep='\t')
        print(f"Count table shape: {count_df.shape}")
        
        # Create a copy for RPKMF calculation
        rpkmf_df = count_df.copy()
        
        # Get contig ID column (either 'Contig' or 'Accession')
        contig_col = 'Contig' if 'Contig' in rpkmf_df.columns else 'Accession'
        
        # Merge with contig lengths
        if contig_col == 'Contig':
            # Strip version numbers for matching (e.g., NC_000898.1 -> NC_000898)
            rpkmf_df['Contig_base'] = rpkmf_df['Contig'].str.replace(r'\.\d+$', '', regex=True)
            merged_lengths = pd.merge(rpkmf_df, contig_lengths, 
                                    left_on='Contig_base', right_on='Accession', how='left')
            # Clean up temporary column
            merged_lengths = merged_lengths.drop('Contig_base', axis=1)
        else:
            # Direct merge on Accession
            merged_lengths = pd.merge(rpkmf_df, contig_lengths, on='Accession', how='left')
        
        print(f"Merged with lengths, missing lengths: {merged_lengths['Length'].isna().sum()}")
        
        # Find count columns (end with '_count')
        count_columns = [col for col in merged_lengths.columns if col.endswith('_count')]
        print(f"Found {len(count_columns)} count columns to convert")
        
        # Convert each count column to RPKMF
        for count_col in count_columns:
            # Extract sample name from column (remove '_count' suffix)
            sample_name = count_col.replace('_count', '')
            
            if sample_name in filtered_reads:
                filtered_read_count = filtered_reads[sample_name]
                
                # RPKMF = (read_count / (contig_length / 1000)) / (filtered_reads / 1e6)
                rpkmf_col = count_col.replace('_count', '_rpkmf')
                
                # Calculate RPKMF
                merged_lengths[rpkmf_col] = (
                    merged_lengths[count_col] / (merged_lengths['Length'] / 1000)
                ) / (filtered_read_count / 1e6)
                
                # Fill NaN values with 0 (for missing lengths or zero counts)
                merged_lengths[rpkmf_col] = merged_lengths[rpkmf_col].fillna(0)
                
            else:
                print(f"Warning: No filtered read count found for sample {sample_name}")
        
        # Drop the Length column that was added during merge (keep original columns)
        if 'Length' in merged_lengths.columns and 'Length' not in count_df.columns:
            merged_lengths = merged_lengths.drop('Length', axis=1)
        
        # Save RPKMF table
        merged_lengths.to_csv(output_path, sep='\t', index=False)
        print(f"Saved RPKMF table: {output_path}")
        print(f"RPKMF table shape: {merged_lengths.shape}")
        
        return merged_lengths
    
    # Convert both tables
    rpkmf_table = convert_table_to_rpkmf(count_table_path, output_rpkmf, "count table")
    rpkmf_with_ref_table = convert_table_to_rpkmf(count_with_ref_path, output_rpkmf_ref, "count with reference table")
    
    print("\n✅ Count to RPKMF conversion completed successfully!")
    
    return rpkmf_table, rpkmf_with_ref_table

def main():
    parser = argparse.ArgumentParser(
        description="Convert CoverM read counts to RPKMF using EsViritu formula"
    )
    
    parser.add_argument(
        "--count-table",
        required=True,
        help="Path to merged count table TSV file"
    )
    
    parser.add_argument(
        "--count-with-ref-table", 
        required=True,
        help="Path to merged count with reference table TSV file"
    )
    
    parser.add_argument(
        "--filtered-reads-dir",
        required=True,
        help="Directory containing filtered read count files"
    )
    
    parser.add_argument(
        "--contig-lengths",
        required=True,
        help="Path to contig lengths TSV file"
    )
    
    parser.add_argument(
        "--output-rpkmf",
        required=True,
        help="Output path for RPKMF table"
    )
    
    parser.add_argument(
        "--output-rpkmf-ref",
        required=True,
        help="Output path for RPKMF with reference table"
    )
    
    args = parser.parse_args()
    
    # Convert count to RPKMF
    convert_count_to_rpkmf(
        count_table_path=args.count_table,
        count_with_ref_path=args.count_with_ref_table,
        filtered_reads_dir=args.filtered_reads_dir,
        contig_lengths_path=args.contig_lengths,
        output_rpkmf=args.output_rpkmf,
        output_rpkmf_ref=args.output_rpkmf_ref
    )

if __name__ == "__main__":
    main()