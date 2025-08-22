#!/bin/bash
# EsViritu Database Setup Script
# Based on: https://github.com/cmmr/EsViritu

echo "Setting up EsViritu database..."

# Create database directory
mkdir -p databases/esviritu_DB
cd databases/esviritu_DB

# Download database (~400 MB)
echo "Downloading EsViritu database (v3.1.1, ~400 MB)..."
wget https://zenodo.org/records/15723755/files/esviritu_db_v3.1.1.tar.gz

# Extract database
echo "Extracting database..."
tar -xvf esviritu_db_v3.1.1.tar.gz

# Clean up
echo "Cleaning up..."
rm esviritu_db_v3.1.1.tar.gz

# Get absolute path
DB_PATH=$(pwd)/v3.1.1

echo "Database setup complete!"
echo "Database location: $DB_PATH"
echo ""
echo "To set the database path in your conda environment:"
echo "conda activate esviritu"
echo "conda env config vars set ESVIRITU_DB=$DB_PATH"
echo ""
echo "Or add this line to your Snakemake shell command:"
echo "export ESVIRITU_DB=$DB_PATH"