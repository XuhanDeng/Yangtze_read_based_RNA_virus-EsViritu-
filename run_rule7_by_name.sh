#!/bin/bash
# Run only Rule 7 by rule name with nohup

# Activate conda environment
conda activate snakemake

# Run only the merge rule by specifying the rule name
nohup snakemake --rule merge_esviritu_results \
    --cores 1 \
    --use-conda \
    --conda-frontend mamba \
    --rerun-incomplete \
    > rule7_byname_$(date +%Y%m%d_%H%M%S).log 2>&1 &

echo "Rule 7 (merge_esviritu_results) submitted to background with nohup"
echo "Check progress with: tail -f rule7_byname_*.log"
echo "Check if running with: jobs"