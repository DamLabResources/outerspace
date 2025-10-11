# UMI Statistics

Unique Molecular Identifiers (UMIs) are short random sequences used to tag and identify unique molecules in sequencing experiments. Statistical analysis of UMI distributions is crucial for:

- Assessing the quality and efficiency of UMI-based experiments
- Identifying potential biases or artifacts in the data
- Understanding the complexity and diversity of the sample
- Validating the effectiveness of UMI error correction

## Configuration-Driven Statistics

The `stats` command uses a configuration-driven approach where each metric is defined separately with its own parameters. The command reads **collapse output** (not count output) and aggregates the data as needed for each statistic.

### Data Aggregation

The stats command works with collapse output and aggregates the data in two ways:

1. **Key occurrence counting** (if only `key_column` is specified): Counts how many times each key value appears. This is useful for counting reads per sequence.

2. **Unique barcode counting** (if both `key_column` and `barcode_column` are specified): Counts how many unique values exist in the barcode column for each key. This is useful for counting unique UMIs per sequence.

### Basic Configuration

Statistics are defined in the configuration file using `[[stats.metrics]]` sections. Each section specifies:

- `method`: The statistical method to use
- `key_column`: Column containing the keys/sequences to analyze
- `barcode_column`: (Optional) Column for counting unique values per key
- `name`: A custom name for the output column
- Additional parameters specific to each method

### Example Configuration

```toml
# Count unique UMIs per corrected sequence
[[stats.metrics]]
method = "simpson_diversity"
key_column = "protospacer_corrected"
barcode_column = "UMI_5prime_UMI_3prime_corrected"
name = "protospacer_simpson"

# Count reads per corrected sequence  
[[stats.metrics]]
method = "shannon_diversity"
key_column = "protospacer_corrected"
name = "protospacer_shannon_reads"

# Calculate error rate by comparing columns
[[stats.metrics]]
method = "error_rate"
original_column = "protospacer"
corrected_column = "protospacer_corrected"
name = "protospacer_error"
```

### Available Methods

- `gini_coefficient`: Gini coefficient calculation
- `shannon_diversity`: Shannon diversity index (supports `base` parameter)
- `simpson_diversity`: Simpson's diversity index
- `umi_recovery_rate`: UMI recovery rate (requires `allowed_list`)
- `umi_efficiency_rate`: UMI efficiency rate (requires `allowed_list`)
- `umi_redundancy`: UMI redundancy metric
- `error_rate`: Error rate from comparing two columns

### Common Parameters

- `key_column`: Column containing the keys to analyze (required for most metrics)
- `barcode_column`: Optional column for grouping/counting
- `allowed_list`: Path to file containing allowed values (required for some metrics)
- `use_corrected`: Whether to use corrected counts (default: true)
- `base`: Logarithm base for Shannon diversity (default: 2.0)
- `original_column`: Original values column (for error_rate)
- `corrected_column`: Corrected values column (for error_rate)

## Key Metrics

Here are important statistical metrics for analyzing UMI distributions:

1. **Gini Coefficient** - Measures the inequality of UMI count distribution. Values range from 0 (perfect equality) to 1 (maximum inequality). Useful for detecting PCR amplification biases.
   - Range: 0 to 1
   - 0: Perfect equality (all UMIs have equal counts)
   - 1: Maximum inequality (one UMI dominates all others)
   - Typical values: 0.3-0.7 for well-behaved experiments

2. **Shannon Diversity Index** - Quantifies the diversity of UMI sequences in a sample, taking into account both richness (number of unique UMIs) and evenness (distribution of counts).
   - Range: 0 to log(N), where N is the number of unique UMIs
   - 0: No diversity (only one UMI present)
   - Higher values: More diverse and even distribution
   - Typical values: 4-8 for well-diversified libraries

3. **Simpson's Diversity Index** - Measures the probability that two randomly selected UMIs are different. More sensitive to dominant UMIs than Shannon's index.
   - Range: 0 to 1
   - 0: No diversity (only one UMI present)
   - 1: Maximum diversity (all UMIs equally abundant)
   - Typical values: 0.7-0.95 for well-diversified libraries


6. **UMI Recovery Rate** - The proportion of expected UMIs that are actually observed after error correction, helping assess the effectiveness of the correction process.
   - Range: 0 to 1
   - 0: No UMIs recovered
   - 1: All expected UMIs recovered
   - Typical values: 0.7-0.9 for good quality data

8. **UMI Error Rate** - The infered sequencing error rate calculated from the merged barcode.
   - Range: 0 to barcode_length
   - 0: No errors

9. **UMI Redundancy** - The average number of reads per unique UMI, useful for assessing the efficiency of the UMI-based deduplication process.
    - Range: 1 to infinity
    - 1: No redundancy (each UMI appears once)
    - Higher values: More redundant sequencing
    - Typical values: 2-10 for well-balanced experiments

These metrics can be used individually or in combination to provide a comprehensive view of UMI-based experiment quality and results. The typical values provided are general guidelines and may vary depending on the specific experimental setup and sequencing technology used.

## Pairwise Library Comparisons

When comparing two UMI libraries (e.g., control vs. treatment, or time point 1 vs. time point 2), several metrics can help quantify their differences and similarities:

1. **Jaccard Similarity** - Measures the overlap between two UMI sets
   - Range: 0 to 1
   - 0: No shared UMIs
   - 1: Identical UMI sets
   - Typical values: 0.1-0.5 for related samples

2. **Bray-Curtis Dissimilarity** - Quantifies the difference in UMI count distributions
   - Range: 0 to 1
   - 0: Identical distributions
   - 1: Completely different distributions
   - Typical values: 0.3-0.7 for related samples

3. **Fold Change** - Ratio of UMI counts between libraries
   - Range: 0 to infinity
   - 1: Equal abundance
   - >1: Enriched in second library
   - <1: Depleted in second library
   - Typical values: 0.1-10 for significant changes

4. **Spearman Correlation** - Measures rank correlation of UMI abundances
   - Range: -1 to 1
   - 1: Perfect positive correlation
   - 0: No correlation
   - -1: Perfect negative correlation
   - Typical values: 0.3-0.8 for related samples

5. **Differential Abundance Score** - Combined measure of fold change and statistical significance
   - Range: -∞ to +∞
   - 0: No significant change
   - Positive: Significant increase
   - Negative: Significant decrease
   - Typical values: -10 to +10 for significant changes

These pairwise metrics are particularly useful for:
- Identifying differentially abundant UMIs
- Assessing the stability of UMI distributions over time
- Detecting systematic biases between libraries
- Validating experimental reproducibility
- Quantifying the magnitude of changes between conditions

When interpreting these metrics, it's important to consider:
- The experimental design and expected changes
- The sequencing depth of each library
- The biological context of the comparison
- The technical variability in the experimental system

## UMI-Level Differential Analysis

When analyzing individual UMIs across experimental groups, several statistical measures can help identify significant changes and their biological relevance:

1. **Log2 Fold Change (LFC)** - Standardized measure of abundance change
   - Range: -∞ to +∞
   - 0: No change
   - Positive: Increased in second group
   - Negative: Decreased in second group
   - Typical values: -3 to +3 for significant changes
   - Interpretation: Each unit represents a 2-fold change

2. **P-value** - Statistical significance of the observed change
   - Range: 0 to 1
   - <0.05: Statistically significant
   - <0.01: Highly significant
   - <0.001: Very highly significant
   - Typical values: <0.05 after multiple testing correction

3. **Effect Size** - Magnitude of the change relative to variability
   - Range: -∞ to +∞
   - 0: No effect
   - |0.2|: Small effect
   - |0.5|: Medium effect
   - |0.8|: Large effect
   - Typical values: |0.5| to |2| for biological changes

4. **Abundance Rank Change** - Change in relative abundance ranking
   - Range: -N to +N (where N is total number of UMIs)
   - 0: No change in rank
   - Positive: Increased in rank
   - Negative: Decreased in rank
   - Typical values: ±100 to ±1000 for significant changes

5. **SEM** - Standard error of the Mean. A measure of reproducibility across replicates
   - Range: 0 to +∞
   - 0: Perfectly consistent
   - +∞: High variability.

These metrics are particularly useful for:
- Identifying differentially abundant UMIs between conditions
- Prioritizing changes for follow-up analysis
- Assessing the reliability of observed changes
- Understanding the magnitude and direction of changes
- Filtering out technical artifacts

When interpreting these metrics, consider:
- The biological context and expected changes
- The number of replicates available
- The sequencing depth in each group
- The overall distribution of changes
- The relationship between statistical and biological significance

Best practices for UMI-level differential analysis:
1. Always use multiple testing correction
2. Consider both statistical and biological significance
3. Validate findings with independent replicates
4. Use appropriate normalization methods
5. Account for batch effects if present
6. Consider the relationship between abundance and variance
7. Use appropriate statistical tests for your experimental design
