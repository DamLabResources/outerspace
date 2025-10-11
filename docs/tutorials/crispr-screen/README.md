# CRISPR Screen Analysis Tutorial

This tutorial demonstrates how to analyze CRISPR screen data using OUTERSPACE's individual commands. We'll process three samples from a CRISPR screen experiment to extract and count barcoded gRNA sequences.

## Overview

We'll analyze three samples:
- `409-4_S1_L002`: Control sample (shuffle)
- `2-G1L9-M1_S9_L001`: First generation library (M1-lib)
- `2-G1L9-M2_S12_L001`: Second generation library (M2-lib)

Each sample has paired-end FASTQ files containing:
- **R1**: 5' UMI + protospacer sequence
- **R2**: 3' UMI

## Prerequisites

1. Clone the OUTERSPACE repository:
```bash
git clone https://github.com/DamLabResources/outerspace.git
cd outerspace
```

2. Install OUTERSPACE:
```bash
pip install -e .
```

3. Navigate to this tutorial directory:
```bash
cd docs/tutorials/crispr-screen
```

## Data Files

This tutorial includes the following data files in the `data/` directory:
- `409-4_S1_L002_R1_001.fastq.gz` / `409-4_S1_L002_R2_001.fastq.gz`
- `2-G1L9-M1_S9_L001_R1_001.fastq.gz` / `2-G1L9-M1_S9_L001_R2_001.fastq.gz`
- `2-G1L9-M2_S12_L001_R1_001.fastq.gz` / `2-G1L9-M2_S12_L001_R2_001.fastq.gz`
- `library_protospacers.txt` (allowed list of expected protospacers)

## Understanding the Configuration File

OUTERSPACE uses TOML configuration files to define patterns and command parameters. Let's break down the `grnaquery.toml` file:

### Global Patterns Section

```toml
[[patterns]]
name = "UMI_5prime"
reg_expr = "(?P<UMI_5prime>.{8})(?:CTTGGCTTTATATATCTTGTGG){s<=4}"
read = "R1"
orientation = "forward"
multiple = "first"
```

This defines a pattern to extract the 5' UMI:
- **`name`**: Identifier for the pattern, becomes a column name in output
- **`reg_expr`**: Regular expression with named capture group `(?P<UMI_5prime>.{8})`
  - `.{8}` captures exactly 8 nucleotides for the UMI
  - `(?:CTTGGCTTTATATATCTTGTGG){s<=4}` matches the constant sequence with up to 4 substitutions
- **`read`**: Which read file to search (R1 or R2)
- **`orientation`**: Search direction (forward or reverse)
- **`multiple`**: How to handle multiple matches (first, last, all)

```toml
[[patterns]]
name = "protospacer"
reg_expr = "(?:TATCTTGTGGAAAGGACGAAACACC){s<=4}(?P<protospacer>.{19,21})(?P<downstreamof_protospacer>GTTTAAGTACTCTGTGCTGGAAACAG){s<=4}"
read = "R1"
orientation = "forward"
multiple = "first"
```

The protospacer pattern extracts the gRNA sequence:
- Matches a constant upstream sequence with fuzzy matching
- Captures the variable protospacer sequence (19-21 nucleotides)
- Also captures a downstream constant region for validation

```toml
[[patterns]]
name = "UMI_3prime"
reg_expr = "(?P<UMI_3prime>.{8})(?:TTCCACACCCTAACTGACACAC){s<=4}"
read = "R2"
orientation = "forward"
multiple = "first"
```

The 3' UMI pattern extracts the second barcode from R2 reads.

### Command-Specific Configurations

```toml
[findseq]
pattern_names = ["UMI_5prime", "protospacer", "UMI_3prime"]
matches_only = true
```

- **`pattern_names`**: Which patterns to use (references the global patterns)
- **`matches_only`**: Only output reads where all patterns are found

```toml
# Iterative collapse steps (recommended approach)
[[collapse.steps]]
name = "umi_correction"
columns = "UMI_5prime,UMI_3prime"
method = "directional"
mismatches = 2

[[collapse.steps]]
name = "protospacer_correction"
columns = "protospacer"
method = "nearest"
allowed_list = "data/library_protospacers.txt"
mismatch_penalty = -1
gap_penalty = -2
match_score = 0
min_score = -3
```

The `[[collapse.steps]]` syntax defines **iterative correction steps**:
- **Step 1 (umi_correction)**: Clusters UMI barcodes using directional method
  - Creates `UMI_5prime_UMI_3prime_corrected` column
- **Step 2 (protospacer_correction)**: Rescues protospacers using nearest-neighbor matching
  - Creates `protospacer_corrected` column
  - Maps near-miss sequences to allowed list entries

This replaces the old workflow where you'd run collapse for UMIs, then use `count --key-rescue` for protospacer correction. Now both corrections happen in a single `collapse` command!

```toml
[count]
barcode_column = 'UMI_5prime_UMI_3prime_corrected'
key_column = 'protospacer_corrected'
```

- **`barcode_column`**: Column containing corrected barcodes to count
- **`key_column`**: Column to group by (now uses the corrected protospacer column)

```toml
[merge]
column = 'UMI_5prime_UMI_3prime_corrected_count'
key_column = 'protospacer_corrected'

[[stats.metrics]]
method = "simpson_diversity"
key_column = "protospacer_corrected"
barcode_column = "UMI_5prime_UMI_3prime_corrected_count"
name = "protospacer_simpson"

[[stats.metrics]]
method = "shannon_diversity"
key_column = "protospacer_corrected"
barcode_column = "UMI_5prime_UMI_3prime_corrected_count"
name = "protospacer_shannon"

[[stats.metrics]]
method = "gini_coefficient"
key_column = "protospacer_corrected"
barcode_column = "UMI_5prime_UMI_3prime_corrected_count"
name = "protospacer_gini"

[[stats.metrics]]
method = "umi_recovery_rate"
key_column = "protospacer_corrected"
barcode_column = "UMI_5prime_UMI_3prime_corrected_count"
allowed_list = "data/library_protospacers.txt"
name = "protospacer_recovery_rate"
```

The `[[stats.metrics]]` sections define individual statistics to calculate:
- Each metric has a `method` (the statistical calculation to perform)
- A `name` for the output column (allowing custom, descriptive names)
- Parameters specific to each metric (e.g., `key_column`, `barcode_column`, `allowed_list`)

This stepwise approach allows you to:
- Calculate multiple different metrics in one run
- Customize which metrics are calculated for your specific needs
- Use different columns or parameters for different metrics
- Get clearly named output columns

## Tutorial Workflow (Configuration File Approach)

This is the **recommended approach** as it ensures consistency and reproducibility across all steps.

### Step 1: Extract Sequences (findseq)

Extract sequences from FASTQ files using the configured patterns:

```bash
# Create output directories
mkdir -p results/findseq results/collapse results/count

# Process control sample (shuffle)
outerspace findseq -c grnaquery.toml \
    -1 data/409-4_S1_L002_R1_001.fastq.gz \
    -2 data/409-4_S1_L002_R2_001.fastq.gz \
    -o results/findseq/shuffle.csv

# Process M1 library sample
outerspace findseq -c grnaquery.toml \
    -1 data/2-G1L9-M1_S9_L001_R1_001.fastq.gz \
    -2 data/2-G1L9-M1_S9_L001_R2_001.fastq.gz \
    -o results/findseq/M1-lib.csv

# Process M2 library sample
outerspace findseq -c grnaquery.toml \
    -1 data/2-G1L9-M2_S12_L001_R1_001.fastq.gz \
    -2 data/2-G1L9-M2_S12_L001_R2_001.fastq.gz \
    -o results/findseq/M2-lib.csv
```

**Expected output**: CSV files containing extracted UMI and protospacer sequences for each sample.

### Step 2: Iterative Collapse (UMI + Protospacer Correction)

Perform multi-stage correction using the iterative steps defined in the config:

```bash
# Perform iterative collapse using steps from config
outerspace collapse -c grnaquery.toml \
    --input-dir results/findseq \
    --output-dir results/collapsed
```

**What this does**:
- **Step 1 (umi_correction)**: Clusters similar UMI barcodes
  - Uses `method = 'directional'` with `mismatches = 2`
  - Creates `UMI_5prime_UMI_3prime_corrected` column
- **Step 2 (protospacer_correction)**: Rescues near-miss protospacers
  - Uses `method = 'nearest'` with allowed list
  - Creates `protospacer_corrected` column
  - Maps sequences like "AAAGGG" to "AAAGGT" if they're similar enough
- Automatically uses temporary directories for intermediate results
- Cleans up temp directories when complete

This replaces the old two-command workflow (collapse + count with key-rescue) with a single streamlined command!

### Step 3: Count Unique Barcodes (count)

Count unique barcodes per corrected protospacer using config settings:

```bash
# Count barcodes for all samples using config
outerspace count -c grnaquery.toml \
    --input-dir results/collapsed \
    --output-dir results/count
```

**What this does**:
- Counts unique `UMI_5prime_UMI_3prime_corrected` barcodes
- Groups by `protospacer_corrected` (already rescued in Step 2)
- No need for `--key-rescue` flag—correction already done!

**Expected output**: CSV files with barcode counts per corrected protospacer for each sample.

### Step 4: Merge Results

Combine all sample results using configuration defaults:

```bash
# Merge in wide format using config
outerspace merge -c grnaquery.toml \
    results/count/shuffle.csv \
    results/count/M1-lib.csv \
    results/count/M2-lib.csv \
    --output-file results/merged_counts_wide.csv \
    --sample-names shuffle M1-lib M2-lib \
    --format wide

# Merge in long format (alternative)
outerspace merge -c grnaquery.toml \
    results/count/shuffle.csv \
    results/count/M1-lib.csv \
    results/count/M2-lib.csv \
    --output-file results/merged_counts_long.csv \
    --sample-names shuffle M1-lib M2-lib \
    --format long
```

### Step 5: Generate Statistics

Calculate summary statistics using the metrics defined in the config:

```bash
# Generate comprehensive statistics using metrics from config
outerspace stats -c grnaquery.toml \
    results/count/shuffle.csv \
    results/count/M1-lib.csv \
    results/count/M2-lib.csv
```

**What this calculates** (based on `[[stats.metrics]]` sections in config):
- **Simpson diversity**: Measures the diversity of protospacer distribution
- **Shannon diversity**: Alternative diversity metric
- **Gini coefficient**: Measures inequality in the distribution
- **UMI redundancy**: Average reads per unique UMI
- **Recovery rate**: Fraction of expected protospacers found
- **Efficiency rate**: Fraction of reads matching allowed protospacers

**Expected output**: CSV file with one row per input file and columns for each metric defined in the config. The output column names match the `name` field from each `[[stats.metrics]]` section.

## Understanding Iterative Collapse Steps

The new `[[collapse.steps]]` feature allows you to chain multiple correction operations in a single command. This is more efficient and cleaner than the old workflow.

### How Protospacer Rescue Works

The second step uses **nearest-neighbor matching** with Needleman-Wunsch alignment:

1. **Exact matches** are kept as-is
2. **Near-miss sequences** are aligned against the allowed list
3. **Alignment scoring** uses match/mismatch/gap penalties
4. **K-mer prescreening** speeds up the search by filtering candidates

**Example**: If the allowed list contains `AAAGGTCTTGCAGCTACGC` and a read has `AAAGTTCTTGCAGCTACGC` (one mismatch), it will be rescued with alignment score = `19*0 + 1*(-1) = -1` (using the tutorial's scoring: match=0, mismatch=-1).

### Tuning Protospacer Rescue Parameters

In your `[[collapse.steps]]` configuration:

```toml
[[collapse.steps]]
name = "protospacer_correction"
columns = "protospacer"
method = "nearest"
allowed_list = "data/library_protospacers.txt"

# Alignment scoring
match_score = 0              # Score for matching bases
mismatch_penalty = -1        # Penalty for mismatches
gap_penalty = -2             # Penalty for insertions/deletions
min_score = -3               # Minimum score to accept (allows ~3 mismatches or 1 indel + 1 mismatch)

# K-mer prescreening (speeds up search)
rescue_kmer_size = 3         # K-mer length for approximate matching
rescue_min_overlap = 8       # Minimum shared k-mers required before alignment
# rescue_exhaustive = true   # Uncomment to disable prescreen (slower but comprehensive)

# Multi-match strategy
rescue_strategy = "random"   # How to choose among tied matches: random, first, last
```

### K-mer Prescreening Strategy

The k-mer prescreen dramatically improves performance:

- **Precomputes** k-mers for all allowed values
- **Filters** candidates that share at least `rescue_min_overlap` k-mers
- **Only aligns** against filtered candidates
- **Trade-off**: Higher `min_overlap` = faster but may miss some rescues

Example: With `rescue_kmer_size = 3` and a sequence `AAAGGT`, the k-mers are: `AAA`, `AAG`, `AGG`, `GGT`. Only allowed sequences sharing ≥8 of these k-mers will be considered for alignment.

## Command-Line Overrides

While configuration files are recommended for reproducible analyses, you can override specific settings via command-line parameters when needed for testing or one-off analyses.

## Single File Processing vs Batch Processing

The `[[collapse.steps]]` feature works with both single file and batch processing modes:

**Batch processing** (recommended for multiple samples):
- Use `--input-dir` and `--output-dir`
- Processes all CSV files in the directory
- Great for processing multiple samples at once

```bash
# Batch mode with iterative steps
outerspace collapse -c grnaquery.toml \
    --input-dir results/findseq/ \
    --output-dir results/collapsed/
```

**Single file processing**:
- Use `--input-file` and `--output-file`
- Also works with `[[collapse.steps]]` for a single sample
- Intermediate results are stored in temporary files and cleaned up automatically

```bash
# Single file mode with iterative steps
outerspace collapse -c grnaquery.toml \
    --input-file results/findseq/shuffle.csv \
    --output-file results/collapsed/shuffle.csv

# Or single file without steps (manual step-by-step)
outerspace collapse \
    --input-file results/findseq/shuffle.csv \
    --output-file results/umi_corrected.csv \
    --columns UMI_5prime,UMI_3prime \
    --mismatches 2 \
    --method directional
```

Both modes benefit from the iterative workflow defined in your config file.

## Understanding the Results

### findseq Output
- Each row represents one sequencing read
- Columns contain extracted sequences: `UMI_5prime`, `protospacer`, `UMI_3prime`
- Only reads matching all patterns are included (due to `matches_only = true`)

### collapsed Output (Iterative Steps)
The iterative collapse creates **two new corrected columns**:

1. **`UMI_5prime_UMI_3prime_corrected`** (from Step 1)
   - Clusters similar UMI barcodes using directional method
   - Reduces sequencing errors in barcodes
   - Example: `AAAACCCC` and `AAAACCCT` → `AAAACCCC`

2. **`protospacer_corrected`** (from Step 2)
   - Maps near-miss sequences to allowed list
   - Uses Needleman-Wunsch alignment
   - Example: `AAAGTTCTTGC...` → `AAAGGTCTTGC...` (rescued 1-mismatch variant)
   - Empty if no rescue possible (score below threshold)

Original columns (`UMI_5prime`, `UMI_3prime`, `protospacer`) are preserved for reference.

### count Output
- Each row represents one unique **corrected** protospacer
- Shows count of unique barcodes per `protospacer_corrected`
- Higher counts indicate more cells with that gRNA
- Only includes successfully rescued/matched protospacers

### merge Output
- **Wide format**: One row per corrected protospacer, one column per sample
- **Long format**: One row per protospacer-sample combination

## Quality Control

Monitor these metrics throughout the analysis:
1. **Extraction efficiency**: Percentage of reads with all patterns found
2. **UMI correction**: Reduction in unique barcodes after Step 1 clustering
3. **Protospacer rescue rate**: Percentage of reads with rescued protospacers (Step 2)
4. **Library representation**: Coverage of expected protospacers in allowed list
5. **Barcode diversity**: Number of unique barcodes per corrected protospacer

The collapse command logs rescue statistics for Step 2, showing how many values were mapped to the allowed list.

## Best Practices

1. **Use `[[collapse.steps]]`** for multi-stage correction pipelines
2. **Use configuration files** for reproducible analyses
3. **Version control your config files** to track parameter changes
4. **Test rescue parameters** with small datasets to optimize `min_score` and k-mer settings
5. **Monitor rescue rates** to ensure appropriate protospacer recovery
6. **Document parameter choices** in your analysis notebooks
7. **Validate patterns** with a subset of data before full processing

## Next Steps

After completing this tutorial, you can:
1. Visualize barcode distributions using plotting tools
2. Perform statistical analysis to identify enriched/depleted gRNAs
3. Calculate Gini coefficients to assess barcode diversity
4. Compare results between different experimental conditions

## Troubleshooting

**Low extraction rates**: Check pattern configurations in `grnaquery.toml`

**High barcode diversity**: Consider adjusting mismatch tolerance in UMI correction step

**Low protospacer rescue rate**: 
- Check alignment scoring parameters (`min_score` may be too stringent)
- Verify allowed list contains expected sequences
- Try `rescue_exhaustive = true` to disable k-mer prescreen
- Lower `rescue_min_overlap` for more candidates

**Too many rescued protospacers** (potential false positives):
- Increase `min_score` threshold
- Increase `rescue_min_overlap` for stricter k-mer filtering
- Review alignment parameters

**Memory issues**: Process samples individually or use row limits for testing

**Steps mode not working**: Ensure you use `--input-dir` and `--output-dir` (not `--input-file`)

**Parameter conflicts**: Command-line arguments always override config file settings


Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
