#!/bin/bash

# CRISPR Screen Analysis Tutorial - Automated Workflow
# This script runs through the complete OUTERSPACE workflow demonstrated in the tutorial
# Uses configuration file approach for consistency and reproducibility

set -e  # Exit on any error

echo "=== OUTERSPACE CRISPR Screen Tutorial ==="
echo "This script demonstrates the configuration file approach for OUTERSPACE workflows."
echo ""

# Create output directories
echo "Creating output directories..."
mkdir -p results/findseq results/collapsed results/count

echo ""
echo "=== Step 1: Extract Sequences (findseq) ==="

# Process control sample (shuffle)
echo "Processing control sample (shuffle)..."
outerspace findseq -c grnaquery.toml \
    -1 data/409-4_S1_L002_R1_001.fastq.gz \
    -2 data/409-4_S1_L002_R2_001.fastq.gz \
    -o results/findseq/shuffle.csv

# Process M1 library sample
echo "Processing M1 library sample..."
outerspace findseq -c grnaquery.toml \
    -1 data/2-G1L9-M1_S9_L001_R1_001.fastq.gz \
    -2 data/2-G1L9-M1_S9_L001_R2_001.fastq.gz \
    -o results/findseq/M1-lib.csv

# Process M2 library sample
echo "Processing M2 library sample..."
outerspace findseq -c grnaquery.toml \
    -1 data/2-G1L9-M2_S12_L001_R1_001.fastq.gz \
    -2 data/2-G1L9-M2_S12_L001_R2_001.fastq.gz \
    -o results/findseq/M2-lib.csv

echo ""
echo "=== Step 2: Iterative Collapse (UMI + Protospacer Correction) ==="

# Perform iterative collapse using steps defined in config
echo "Running iterative collapse (2 steps):"
echo "  Step 1: UMI correction using directional clustering"
echo "  Step 2: Protospacer correction using nearest-neighbor matching"
outerspace collapse -c grnaquery.toml \
    --input-dir results/findseq \
    --output-dir results/collapsed

echo ""
echo "=== Step 3: Count Unique Barcodes (count) ==="

# Count barcodes for all samples using config
# Note: Protospacer correction was already done in Step 2
echo "Counting barcodes for all samples using configuration file..."
outerspace count -c grnaquery.toml \
    --input-dir results/collapsed \
    --output-dir results/count

echo ""
echo "=== Step 4: Merge Results ==="

# Merge in wide format using config
echo "Merging results in wide format using configuration file..."
outerspace merge -c grnaquery.toml \
    results/count/shuffle.csv \
    results/count/M1-lib.csv \
    results/count/M2-lib.csv \
    --output-file results/merged_counts_wide.csv \
    --sample-names shuffle M1-lib M2-lib \
    --format wide

# Merge in long format using config
echo "Merging results in long format using configuration file..."
outerspace merge -c grnaquery.toml \
    results/count/shuffle.csv \
    results/count/M1-lib.csv \
    results/count/M2-lib.csv \
    --output-file results/merged_counts_long.csv \
    --sample-names shuffle M1-lib M2-lib \
    --format long

echo ""
echo "=== Step 5: Generate Statistics ==="

# Generate comprehensive statistics using config
echo "Generating statistics for all samples using configuration file..."
outerspace stats -c grnaquery.toml \
    results/collapsed/shuffle.csv \
    results/collapsed/M1-lib.csv \
    results/collapsed/M2-lib.csv >  results/statistics.csv

echo ""
echo "=== Tutorial Complete! ==="
echo ""
echo "This tutorial demonstrated the iterative collapse workflow with:"
echo "  ✓ Multi-step correction pipeline defined in config file"
echo "  ✓ Automatic UMI correction followed by protospacer rescue"
echo "  ✓ Consistent parameter usage across all commands"
echo "  ✓ Easy reproducibility with version-controlled config files"
echo ""
echo "Results have been saved to the 'results/' directory:"
echo "  - results/findseq/: Extracted sequences"
echo "  - results/collapsed/: Iteratively corrected files (UMI + protospacer)"
echo "    ├─ UMI_5prime_UMI_3prime_corrected: Clustered UMI barcodes"
echo "    └─ protospacer_corrected: Rescued protospacer sequences"
echo "  - results/count/: Barcode counts per corrected protospacer"
echo "  - results/merged_*.csv: Merged results across samples"
echo ""
echo "Key takeaways:"
echo "  • [[collapse.steps]] enables iterative multi-stage correction"
echo "  • Combines UMI clustering + nearest-neighbor rescue in one command"
echo "  • Configuration files ensure consistency and reproducibility"
echo "  • Version control your config files for better analysis tracking"
echo ""
echo "You can now explore the results and run additional analyses." 


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.