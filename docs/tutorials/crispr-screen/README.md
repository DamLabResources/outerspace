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
[collapse]
columns = 'UMI_5prime,UMI_3prime'
mismatches = 2
method = 'directional'
```

- **`columns`**: Which columns to combine for barcode correction
- **`mismatches`**: Maximum edit distance for clustering similar barcodes
- **`method`**: UMI-tools clustering algorithm

```toml
[count]
barcode_column = 'UMI_5prime_UMI_3prime_corrected'
key_column = 'protospacer'
```

- **`barcode_column`**: Column containing corrected barcodes to count
- **`key_column`**: Column to group by (protospacer sequences)

```toml
[merge]
column = 'UMI_5prime_UMI_3prime_corrected_count'
key_column = 'protospacer'

[stats]
key_column = 'protospacer'
count_column = 'UMI_5prime_UMI_3prime_corrected_count'
```

These sections define default parameters for the merge and stats commands.

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

### Step 2: Correct Barcodes (collapse)

Correct sequencing errors using the configuration file settings:

```bash
# Correct barcodes for all samples using config
outerspace collapse -c grnaquery.toml \
    --input-dir results/findseq \
    --output-dir results/collapse
```

**What this does**:
- Uses `columns = 'UMI_5prime,UMI_3prime'` from config
- Applies `mismatches = 2` and `method = 'directional'` from config
- Creates corrected barcode column: `UMI_5prime_UMI_3prime_corrected`

### Step 3: Count Unique Barcodes (count)

Count unique barcodes per protospacer using config settings:

```bash
# Count barcodes for all samples using config
outerspace count -c grnaquery.toml \
    --input-dir results/collapse \
    --output-dir results/count
```

**Expected output**: CSV files with barcode counts per protospacer for each sample.

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

Calculate summary statistics using config settings:

```bash
# Generate comprehensive statistics using config
outerspace stats -c grnaquery.toml \
    results/count/shuffle.csv \
    results/count/M1-lib.csv \
    results/count/M2-lib.csv
```

## Advanced Workflow: Using Allowed Lists

For more stringent analysis, you can filter results to only include expected protospacers:

### Step 3b: Count with Allowed List and Key Rescue

```bash
# Count barcodes using only library protospacers and rescue near-miss keys
outerspace count -c grnaquery.toml \
    --input-dir results/collapse \
    --output-dir results/count_filtered \
    --allowed-list data/library_protospacers.txt \
    --key-rescue --key-min-score 17 --key-match-score 1 \
    --key-mismatch-penalty -1 --key-gap-penalty -3
```

#### Prescreening strategy for key rescue (k-mer approximate matching)

When `--key-rescue` is enabled, OUTERSPACE accelerates nearest-key lookup with a k-mer prescreen:

- It precomputes sliding k-mers for every allowed key
- For each query key, it computes its k-mers and only aligns against allowed keys that share at least a minimum number of k-mers
- This radically reduces expensive alignments without changing final results under typical settings

You can tune the prescreen with these flags:

- `--key-rescue-kmer-size`: k-mer length used for prescreening (default: 3)
- `--key-rescue-min-overlap`: minimum shared k-mers required to attempt alignment (default: 1)
- `--key-rescue-exhaustive`: disable prescreen and align against all allowed keys (slower)

Examples:

```bash
# Stricter prescreen (fewer alignments)
outerspace count -c grnaquery.toml \
  --input-dir results/collapse \
  --output-dir results/count_filtered \
  --allowed-list data/library_protospacers.txt \
  --key-rescue \
  --key-rescue-kmer-size 4 \
  --key-rescue-min-overlap 2
```

```bash
# Exhaustive search (no prescreen)
outerspace count -c grnaquery.toml \
  --input-dir results/collapse \
  --output-dir results/count_filtered \
  --allowed-list data/library_protospacers.txt \
  --key-rescue \
  --key-rescue-exhaustive
```

In TOML configs, you can set:

```toml
[count]
key_rescue = true
key_rescue_kmer_size = 3     # same as --key-rescue-kmer-size
key_rescue_min_overlap = 8   # same as --key-rescue-min-overlap
```

Command-line flags override TOML values when both are provided. Use exhaustive mode from the CLI if you want to temporarily disable the prescreen for testing.

**Notes**:
- The `--allowed-list` parameter enables filtering to expected protospacers.
- You can enable key rescue directly in the TOML (recommended for increasing the number of useful reads):

```toml
[count]
allowed_list = 'data/library_protospacers.txt'
key_rescue = true
key_mismatch_penalty = -1
key_gap_penalty = -3
key_match_score = 1
key_min_score = 17
```

Command-line flags will override TOML values when provided.

### Step 4b: Merge Filtered Results

```bash
mkdir -p results/count_filtered

# Merge filtered results
outerspace merge -c grnaquery.toml \
    results/count_filtered/shuffle.csv \
    results/count_filtered/M1-lib.csv \
    results/count_filtered/M2-lib.csv \
    --output-file results/merged_filtered_counts.csv \
    --sample-names shuffle M1-lib M2-lib \
    --format wide
```

## Command-Line Overrides

While configuration files are recommended for reproducible analyses, you can override specific settings via command-line parameters when needed for testing or one-off analyses.

## Single File Processing

You can also process files individually instead of using directories:

```bash
# Process a single sample through the entire workflow
outerspace findseq -c grnaquery.toml \
    -1 data/409-4_S1_L002_R1_001.fastq.gz \
    -2 data/409-4_S1_L002_R2_001.fastq.gz \
    -o results/single_shuffle.csv

outerspace collapse -c grnaquery.toml \
    --input-file results/single_shuffle.csv \
    --output-file results/single_shuffle_collapsed.csv

outerspace count -c grnaquery.toml \
    --input-file results/single_shuffle_collapsed.csv \
    --output-file results/single_shuffle_counts.csv
```

## Understanding the Results

### findseq Output
- Each row represents one sequencing read
- Columns contain extracted sequences (UMI_5prime, protospacer, UMI_3prime)
- Only reads matching all patterns are included (due to `matches_only = true`)

### collapse Output
- Adds corrected barcode columns (e.g., `UMI_5prime_UMI_3prime_corrected`)
- Groups similar barcodes to reduce sequencing errors
- Preserves original barcode information

### count Output
- Each row represents one unique protospacer
- Shows count of unique barcodes per protospacer
- Higher counts indicate more cells with that gRNA

### merge Output
- **Wide format**: One row per protospacer, one column per sample
- **Long format**: One row per protospacer-sample combination

## Quality Control

Monitor these metrics throughout the analysis:
1. **Extraction efficiency**: Percentage of reads with all patterns found
2. **Barcode diversity**: Number of unique barcodes per protospacer
3. **Library representation**: Coverage of expected protospacers
4. **Barcode correction**: Reduction in unique barcodes after clustering

## Best Practices

1. **Use configuration files** for reproducible analyses
2. **Version control your config files** to track parameter changes
3. **Test with small datasets** before processing large files
4. **Document parameter choices** in your analysis notebooks
5. **Validate patterns** with a subset of data before full processing

## Next Steps

After completing this tutorial, you can:
1. Visualize barcode distributions using plotting tools
2. Perform statistical analysis to identify enriched/depleted gRNAs
3. Calculate Gini coefficients to assess barcode diversity
4. Compare results between different experimental conditions

## Troubleshooting

**Low extraction rates**: Check pattern configurations in `grnaquery.toml`
**High barcode diversity**: Consider adjusting mismatch tolerance in collapse step
**Missing protospacers**: Verify allowed list contains expected sequences
**Memory issues**: Process samples individually or use row limits for testing
**Parameter conflicts**: Command-line arguments always override config file settings


Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
