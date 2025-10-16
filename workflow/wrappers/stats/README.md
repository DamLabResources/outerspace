# Stats Wrapper

This wrapper provides a Snakemake interface for the outerspace stats command, which calculates statistics from UMI count data using config-defined metrics.

## Input/Output

### Input
- `csv`: CSV file(s) containing collapsed UMI data to analyze (can be single file or list of files)
- `toml`: TOML configuration file defining the metrics to calculate (see `[[stats.metrics]]` sections)

### Output
- `output`: CSV file containing calculated statistics for all input files

### Threads
- `threads`: Number of threads for parallel file processing (optional, default: 1)

## Example

```python
# Single file
rule stats_single:
    input:
        csv = "collapse/sample1.csv",
        toml = "config.toml"
    output:
        "stats/sample1_stats.csv"
    threads: 1
    wrapper:
        "wrappers/stats"

# Multiple files with threading
rule stats_multiple:
    input:
        csv = ["collapse/sample1.csv", "collapse/sample2.csv", "collapse/sample3.csv"],
        toml = "config.toml"
    output:
        "stats/all_stats.csv"
    threads: 4
    wrapper:
        "wrappers/stats"
```

## Configuration

Statistics are defined in the TOML config file using `[[stats.metrics]]` sections. Each metric specifies:
- `method`: The statistical method to use (e.g., "gini_coefficient", "shannon_diversity", "hill_number")
- `name`: Output column name for this metric
- Additional metric-specific parameters (e.g., `key_column`, `barcode_column`, `q` for Hill numbers)

Example TOML config:
```toml
[[stats.metrics]]
method = "gini_coefficient"
name = "gini"
key_column = "protospacer_corrected"
barcode_column = "UMI_5prime_UMI_3prime_corrected"

[[stats.metrics]]
method = "hill_number"
name = "hill_shannon"
key_column = "protospacer_corrected"
barcode_column = "UMI_5prime_UMI_3prime_corrected"
q = "shannon"
```

## Notes

- The TOML config file must contain `[[stats.metrics]]` sections defining which metrics to calculate
- Available metrics: gini_coefficient, shannon_diversity, simpson_diversity, hill_number, umi_recovery_rate, umi_efficiency_rate, error_rate
- When multiple files are provided, statistics for each file are included as separate rows in the output CSV
- The output CSV includes a 'filename' column to identify which file each row corresponds to
- Using `threads > 1` enables parallel processing of multiple input files for faster computation 


### Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.