#!/bin/bash
# Run Rule 5 and Rule 7 by rule names with nohup

# Activate conda environment
conda activate snakemake

# Run both esviritu identification and merge rules by specifying rule names
nohup snakemake --rule esviritu_identification --rule merge_esviritu_results \
    --cores 8 \
    --use-conda \
--conda-frontend mamba \
    --rerun-incomplete \
    > rule5_and_rule7_$(date +%Y%m%d_%H%M%S).log 2>&1 &

echo "Rule 5 (esviritu_identification) and Rule 7 (merge_esviritu_results) submitted to background with nohup"
echo "Check progress with: tail -f rule5_and_rule7_*.log"
echo "Check if running with: jobs"