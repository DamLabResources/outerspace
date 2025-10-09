# Iterative Collapse with Steps

The `collapse` command now supports defining iterative correction steps in your TOML configuration file. This allows you to chain multiple collapse operations together, where the output of one step becomes the input for the next.

## Configuration Format

Define steps using the `[[collapse.steps]]` array in your TOML config:

```toml
[[collapse.steps]]
name = "step_name"
columns = "column1,column2"
method = "directional"
mismatches = 2
# ... other parameters

[[collapse.steps]]
name = "another_step"
columns = "column3"
method = "nearest"
allowed_list = "allowed_values.txt"
# ... other parameters
```

## Usage

When steps are defined in the config, use the `collapse` command with either single file or batch processing:

**Batch processing** (multiple files):
```bash
outerspace collapse \
  --config my_config.toml \
  --input-dir input_data/ \
  --output-dir final_output/
```

**Single file processing**:
```bash
outerspace collapse \
  --config my_config.toml \
  --input-file input_data/sample.csv \
  --output-file final_output/sample.csv
```

The command will:
1. Execute each step in sequence
2. Use temporary directories (batch) or temporary files (single file) for intermediate results
3. Clean up temp paths automatically
4. Place final results in the specified output location

## Example: UMI + Protospacer Correction

```toml
# Step 1: Correct UMIs using directional clustering
[[collapse.steps]]
name = "umi_correction"
columns = "UMI_5prime,UMI_3prime"
method = "directional"
mismatches = 2

# Step 2: Correct protospacers using nearest-neighbor matching
[[collapse.steps]]
name = "protospacer_correction"
columns = "protospacer"
method = "nearest"
allowed_list = "library_protospacers.txt"
mismatch_penalty = -1
gap_penalty = -3
match_score = 1
min_score = 0
```

This configuration:
- First creates `UMI_5prime_UMI_3prime_corrected` column
- Then creates `protospacer_corrected` column using the allowed list

## Available Step Parameters

### Common Parameters
- `name`: Descriptive name for the step (optional)
- `columns`: Comma-separated list of columns to correct
- `method`: Correction method (`cluster`, `adjacency`, `directional`, `allowed`, `nearest`)
- `sep`: CSV separator (default: `,`)

### UMI Clustering Methods (`cluster`, `adjacency`, `directional`)
- `mismatches`: Number of mismatches allowed (default: 2)

### Allowed List Methods (`allowed`, `nearest`)
- `allowed_list`: Path to file containing allowed values (required)

### Nearest Method Only
- `mismatch_penalty`: Score penalty for mismatches (default: -1)
- `gap_penalty`: Score penalty for gaps/indels (default: -3)
- `match_score`: Score for matches (default: 1)
- `min_score`: Minimum score to accept a match (default: 0)
- `rescue_exhaustive`: Disable k-mer prescreen for exhaustive search (default: false)
- `rescue_kmer_size`: K-mer size for prescreening (default: 3)
- `rescue_min_overlap`: Minimum k-mer overlap required (default: 1)
- `rescue_strategy`: Strategy for multiple matches: `random`, `first`, `last` (default: `random`)

## Relative Paths

The `allowed_list` path can be relative to the config file location:

```toml
# Config at: /path/to/configs/analysis.toml
[[collapse.steps]]
allowed_list = "../data/allowed_keys.txt"  # Resolves to /path/to/data/allowed_keys.txt
```

## Benefits

1. **Cleaner Workflows**: Define entire pipeline in one config file
2. **Automatic Cleanup**: Intermediate results are automatically removed
3. **Better Organization**: Each step is clearly documented
4. **Flexible**: Mix different correction methods as needed

