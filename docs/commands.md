# Commands

Below are the commands available in `outerspace`:

### `outerspace findseq`
Extracts sequences from various file formats based on configuration patterns. Features:
- Support for FASTQ, FASTA, SAM, and BAM files
- Single file or paired-end processing
- Global pattern configuration system
- Auto-detection of file formats
- Multi-threaded processing support
- Progress tracking and error handling

```bash
usage: outerspace findseq [-h] [-1 READ1_FILENAME] [-2 READ2_FILENAME] 
                          [-o OUTPUT_FILENAME] [--region REGION] 
                          [--fetch {mapped,unmapped,all}] [--long-format] 
                          [--matches-only] [--threads THREADS] [--skip-unmapped] 
                          [--max-reads MAX_READS] [--config CONFIG] 
                          [--progress-bar] [--log-file LOG_FILE] [--log-level LOG_LEVEL]

Extract sequences from files based on configuration patterns

options:
  -h, --help            show this help message and exit
  -1 READ1_FILENAME, --read1_filename READ1_FILENAME
                        Input file for read 1 (FASTQ, FASTA, SAM, BAM) or single read file
  -2 READ2_FILENAME, --read2_filename READ2_FILENAME
                        Input file for read 2 (FASTQ, FASTA, SAM, BAM) for paired reads
  -o OUTPUT_FILENAME, --output_filename OUTPUT_FILENAME
                        Output CSV file name
  --region REGION       SAM/BAM region specification (e.g., "chr1:1-1000")
  --fetch {mapped,unmapped,all}
                        SAM/BAM fetch mode (mapped, unmapped, all)
  --long-format         Output in long format (one row per pattern match instead of one row per read)
  --matches-only        Only output reads that have at least one pattern match
  --threads THREADS     Number of threads for parallel processing (default: auto-detect)
  --skip-unmapped       Skip unmapped reads in SAM/BAM files (based on SAM flag 0x4)
  --max-reads MAX_READS
                        Maximum number of reads to process (default: all)

Common Arguments:
  --config CONFIG, -c CONFIG
                        Configuration file
  --progress-bar, -p    Enable progress bar
  --log-file LOG_FILE   Log file
  --log-level LOG_LEVEL
                        Log level (default: WARNING)
```

### `outerspace collapse`
Corrects barcodes in CSV files using UMI-tools clustering. Features:
- Supports multiple barcode columns
- Configurable mismatch tolerance
- Multiple clustering methods (cluster, adjacency, directional)
- Row limiting for testing
- Detailed metrics output

```bash
usage: outerspace collapse [-h] (--input-file INPUT_FILE | --input-dir INPUT_DIR)
                           (--output-file OUTPUT_FILE | --output-dir OUTPUT_DIR)
                           [--columns COLUMNS] [--mismatches MISMATCHES] [--sep SEP]
                           [--row-limit ROW_LIMIT] [--method {cluster,adjacency,directional}]
                           [--metrics METRICS] [--config CONFIG] [--progress-bar]
                           [--log-file LOG_FILE] [--log-level LOG_LEVEL]

Correct barcodes in CSV files using UMI-tools clustering

options:
  -h, --help            show this help message and exit
  --input-file INPUT_FILE
                        Input CSV file to process
  --input-dir INPUT_DIR
                        Input directory containing CSV files to process
  --output-file OUTPUT_FILE
                        Output CSV file for corrected barcodes
  --output-dir OUTPUT_DIR
                        Output directory for corrected CSV files
  --columns COLUMNS     Column(s) containing barcodes to correct. Can be a single column or comma-separated list
  --mismatches MISMATCHES
                        Number of mismatches allowed for clustering (default: 2)
  --sep SEP             CSV separator (default: ,)
  --row-limit ROW_LIMIT
                        Process only the first N rows (for testing)
  --method {cluster,adjacency,directional}
                        Clustering method to use (default: directional)
  --metrics METRICS     Output YAML file for metrics

Common Arguments:
  --config CONFIG, -c CONFIG
                        Configuration file
  --progress-bar, -p    Enable progress bar
  --log-file LOG_FILE   Log file
  --log-level LOG_LEVEL
                        Log level (default: WARNING)
```

### `outerspace count`
Counts unique barcodes per key value in CSV files. Features:
- Barcode and key column specification
- `--allowed-list` filtering
- Optional key rescue with edit-tolerant matching to allowed keys
- Downsampling capability
- Detailed output with barcode lists
- Gini coefficient calculation for both barcodes and keys

```bash
usage: outerspace count [-h] (--input-file INPUT_FILE | --input-dir INPUT_DIR)
                        (--output-file OUTPUT_FILE | --output-dir OUTPUT_DIR)
                        [--barcode-column BARCODE_COLUMN] [--key-column KEY_COLUMN]
                        [--sep SEP] [--row-limit ROW_LIMIT] [--allowed-list ALLOWED_LIST]
                        [--key-rescue] [--key-min-score KEY_MIN_SCORE]
                        [--key-match-score KEY_MATCH_SCORE]
                        [--key-mismatch-penalty KEY_MISMATCH_PENALTY]
                        [--key-gap-penalty KEY_GAP_PENALTY]
                        [--key-rescue-strategy KEY_RESCUE_STRATEGY]
                        [--detailed] [--downsample DOWNSAMPLE] [--random-seed RANDOM_SEED]
                        [--config CONFIG] [--progress-bar] [--log-file LOG_FILE] 
                        [--log-level LOG_LEVEL]

Count unique barcodes per key value in CSV files

options:
  -h, --help            show this help message and exit
  --input-file INPUT_FILE
                        Input CSV file to process
  --input-dir INPUT_DIR
                        Input directory containing CSV files to process
  --output-file OUTPUT_FILE
                        Output CSV file for barcode counts
  --output-dir OUTPUT_DIR
                        Output directory for barcode counts
  --barcode-column BARCODE_COLUMN
                        Column containing barcodes
  --key-column KEY_COLUMN
                        Column to group by
  --sep SEP             CSV separator (default: ,)
  --row-limit ROW_LIMIT
                        Process only the first N rows (for testing)
  --allowed-list ALLOWED_LIST
                        Text file containing allowed keys (one per line)
  --key-rescue          If --allowed-list is provided, rescue near-miss keys by
                        aligning to the closest allowed key
  --key-min-score KEY_MIN_SCORE
                        Minimum alignment score required to accept a rescued key (default: 0)
  --key-match-score KEY_MATCH_SCORE
                        Score for character matches when aligning keys (default: 1)
  --key-mismatch-penalty KEY_MISMATCH_PENALTY
                        Penalty for mismatches when aligning keys (default: -1)
  --key-gap-penalty KEY_GAP_PENALTY
                        Penalty for gaps/indels when aligning keys (default: -3)
  --key-rescue-strategy KEY_RESCUE_STRATEGY
                        Strategy to choose among multiple equally good rescued keys
                        (choices: random, first, last, all; default: random)
  --detailed            Include barcode lists in output
  --downsample DOWNSAMPLE
                        Randomly sample reads with probability between 0 and 1
  --random-seed RANDOM_SEED
                        Random seed for downsampling

Common Arguments:
  --config CONFIG, -c CONFIG
                        Configuration file
  --progress-bar, -p    Enable progress bar
  --log-file LOG_FILE   Log file
  --log-level LOG_LEVEL
                        Log level (default: WARNING)
```

### `outerspace merge`
Merges multiple UMI count files into a single file. Features:
- Wide and long output formats
- Optional UMI clustering with various methods
- Comprehensive metrics reporting
- Sample name customization

```bash
usage: outerspace merge [-h] [--output-file OUTPUT_FILE] [--key-column KEY_COLUMN]
                        [--count-column COUNT_COLUMN] [--sample-names SAMPLE_NAMES [SAMPLE_NAMES ...]]
                        [--sep SEP] [--format {wide,long}] [--mismatches MISMATCHES]
                        [--method {cluster,adjacency,directional}] [--metrics METRICS]
                        [--config CONFIG] [--progress-bar] [--log-file LOG_FILE]
                        [--log-level LOG_LEVEL] files [files ...]

Merge multiple UMI count files into a single file

positional arguments:
  files                 Input CSV files to merge

options:
  -h, --help            show this help message and exit
  --output-file OUTPUT_FILE
                        Output CSV file for merged counts
  --key-column KEY_COLUMN
                        Column containing UMIs
  --count-column COUNT_COLUMN
                        Column containing counts (if not provided, assumes count=1)
  --sample-names SAMPLE_NAMES [SAMPLE_NAMES ...]
                        Optional list of sample names (must match number of input files)
  --sep SEP             CSV separator (default: ,)
  --format {wide,long}  Output format: wide (samples as columns) or long (sample,umi,count columns) (default: wide)
  --mismatches MISMATCHES
                        Number of mismatches allowed for clustering (default: 0)
  --method {cluster,adjacency,directional}
                        Clustering method to use (default: directional)
  --metrics METRICS     Output YAML file for metrics

Common Arguments:
  --config CONFIG, -c CONFIG
                        Configuration file
  --progress-bar, -p    Enable progress bar
  --log-file LOG_FILE   Log file
  --log-level LOG_LEVEL
                        Log level (default: WARNING)
```

### `outerspace stats`
Calculates comprehensive statistics from UMI count data. Features:
- Multiple diversity metrics (Gini coefficient, Shannon diversity, Simpson diversity)
- Efficiency measures (recovery rate, efficiency rate, error rate)
- Redundancy analysis
- Support for pre-counted values
- Scale factor application for normalized counts
- `--allowed-list` filtering

```bash
usage: outerspace stats [-h] [--key-column KEY_COLUMN] [--count-column COUNT_COLUMN]
                        [--scale SCALE] [--sep SEP] [--allowed-list ALLOWED_LIST]
                        [--config CONFIG] [--progress-bar] [--log-file LOG_FILE]
                        [--log-level LOG_LEVEL] input_files [input_files ...]

Calculate all single-library statistics from counts in a CSV column

positional arguments:
  input_files           Input CSV file(s) to process (supports glob patterns)

options:
  -h, --help            show this help message and exit
  --key-column KEY_COLUMN
                        Column containing keys
  --count-column COUNT_COLUMN
                        Column containing pre-counted values
  --scale SCALE         Scale factor for normalized values (e.g., if normalized to mean=1)
  --sep SEP             CSV separator (default: ,)
  --allowed-list ALLOWED_LIST
                        Text file containing allowed values (one per line)

Common Arguments:
  --config CONFIG, -c CONFIG
                        Configuration file
  --progress-bar, -p    Enable progress bar
  --log-file LOG_FILE   Log file
  --log-level LOG_LEVEL
                        Log level (default: WARNING)
```

### `outerspace pipeline`
Runs complete OUTERSPACE workflows using Snakemake. Features:
- Integrated workflow execution
- Configuration management
- Comprehensive error handling
- Custom Snakemake argument support

```bash
usage: outerspace pipeline [-h] [--snakemake-args SNAKEMAKE_ARGS] 
                           config_file snakemake_config

Run the complete OUTERSPACE pipeline using Snakemake

positional arguments:
  config_file           TOML configuration file with search patterns
  snakemake_config      YAML configuration file for Snakemake workflow

options:
  -h, --help            show this help message and exit
  --snakemake-args SNAKEMAKE_ARGS
                        Additional arguments to pass to Snakemake (e.g. --snakemake-args="--dry-run --cores 4")
```

### `outerspace visualize`
Creates visualizations of barcode counts from CSV files. Features:
- Histogram generation
- Customizable plot parameters
- Multiple output formats
- Log scale support

```bash
usage: outerspace visualize [-h] [--sep SEP] [--bins BINS] [--title-prefix TITLE_PREFIX]
                            [--xlabel XLABEL] [--ylabel YLABEL] [--log-scale]
                            [--format FORMAT] [--config CONFIG] [--progress-bar]
                            [--log-file LOG_FILE] [--log-level LOG_LEVEL]
                            input_dir output_dir

Visualize barcode counts from CSV files

positional arguments:
  input_dir             Input directory containing CSV files with barcode counts
  output_dir            Output directory for visualization plots

options:
  -h, --help            show this help message and exit
  --sep SEP             CSV separator (default: ,)
  --bins BINS           Number of histogram bins (default: 50)
  --title-prefix TITLE_PREFIX
                        Prefix for plot titles (default: filename)
  --xlabel XLABEL       X-axis label (default: Number of Unique Barcodes)
  --ylabel YLABEL       Y-axis label (default: Count)
  --log-scale           Use log scale for y-axis
  --format FORMAT       Output image format (default: png)

Common Arguments:
  --config CONFIG, -c CONFIG
                        Configuration file
  --progress-bar, -p    Enable progress bar
  --log-file LOG_FILE   Log file
  --log-level LOG_LEVEL
                        Log level (default: WARNING)
```

Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.