#!/bin/bash
# Test database download rule locally

conda activate snakemake

# Clean any existing attempts
rm -rf databases/esviritu_DB log/0_database

# Test the download rule
snakemake download_esviritu_database --cores 1 --use-conda --dry-run

echo "Dry run completed. If no errors, run without --dry-run"