# Snakemake Integration

## Overview

Snakemake is a workflow management system that helps create reproducible and scalable data analyses. It allows you to define workflows using a Python-based language and automatically handles job scheduling, parallelization, and dependency management.

## Pipeline Command

The `outerspace pipeline` command provides a thin wrapper around Snakemake to run the OUTERSPACE analysis workflow. It takes two main inputs:
1. A TOML configuration file containing search patterns and analysis parameters
2. A YAML configuration file for Snakemake workflow settings
3. And optionally a set of arguments to pass to the Snakemake command

The command automatically:
- Loads and merges configurations
- Sets up the Snakemake environment
- Handles additional Snakemake arguments (including executor selection)
- Provides proper error handling and logging

### Basic Usage

```bash
# Local execution with 4 cores
outerspace pipeline config.toml snakemake_config.yaml --snakemake-args="--cores 4"

# Dry run to test the workflow
outerspace pipeline config.toml snakemake_config.yaml --snakemake-args="--dry-run"
```

### Using SLURM Executor

To run the pipeline on a SLURM cluster, you need to:

1. Install the SLURM executor plugin:
   ```bash
   pip install snakemake-executor-plugin-slurm
   ```

2. Run with the SLURM executor:
   ```bash
   outerspace pipeline config.toml snakemake_config.yaml \
       --snakemake-args="--executor slurm --jobs 100"
   ```

3. (Optional) Specify SLURM-specific resources:
   ```bash
   outerspace pipeline config.toml snakemake_config.yaml \
       --snakemake-args="--executor slurm --jobs 100 \
       --default-resources slurm_account=myaccount slurm_partition=compute"
   ```

**Important Notes:**
- The pipeline automatically creates a `.snakemake/log` directory for SLURM logs
- SLURM job logs will be written to `.snakemake/log/` by default
- Make sure you have sufficient permissions in the working directory
- For large workflows, consider using `--latency-wait 60` to account for filesystem delays

### Using Profiles

Snakemake profiles allow you to store commonly-used settings in a configuration file, making it easier to run workflows consistently.

**Creating a Profile:**

Create a directory for your profile with a `config.yaml` or `config.v8+.yaml` file:

```bash
mkdir -p profiles/slurm
```

Create `profiles/slurm/config.v8+.yaml`:
```yaml
# Executor settings
executor: slurm
jobs: 100

# SLURM-specific settings (passed to the executor plugin)
slurm_partition: compute
slurm_account: myproject
slurm_qos: normal

# Resource defaults (applied to all jobs unless overridden)
# IMPORTANT: These are crucial for SLURM execution
# The slurm_account is often REQUIRED by cluster configurations
# You can use either list format (shown here) or dict format (see below)
default-resources:
  - slurm_account=myproject      # REQUIRED for most SLURM clusters
  - slurm_partition=compute      # Which partition to submit to
  - mem_mb=4000                  # Default memory per job
  - runtime=120                  # Default runtime in minutes
```

**Alternative dict format** (both formats are supported):
```yaml
default-resources:
  slurm_account: myproject
  slurm_partition: compute
  mem_mb: 4000
  runtime: 120
```

**Using a Profile:**

```bash
outerspace pipeline config.toml snakemake_config.yaml \
    --snakemake-args="--profile profiles/slurm"
```

**Profile Priority:**

Settings are applied in this order (later overrides earlier):
1. Profile settings
2. Command-line arguments

For example, this will use the profile but override jobs:
```bash
outerspace pipeline config.toml snakemake_config.yaml \
    --snakemake-args="--profile profiles/slurm --jobs 200"
```

### Using Other Executors

Snakemake v9 supports various executor plugins. Common options include:

- **Local**: `--executor local --cores 8` (default, runs on local machine)
- **SLURM**: `--executor slurm --jobs 100` (SLURM cluster)
- **LSF**: `--executor lsf --jobs 100` (IBM LSF cluster)
- **Grid Engine**: `--executor cluster-generic --jobs 100` (SGE/UGE)
- **Google Cloud**: `--executor googlebatch --jobs 100`

For executor-specific options, refer to the respective executor plugin documentation.

Note: Executor plugins must be installed separately via pip (e.g., `pip install snakemake-executor-plugin-slurm`).

## Troubleshooting

**Error: SLURM jobs not submitting**

Check that:
- You have access to the SLURM partition you're requesting
- Your account has the necessary permissions
- Resource requests are within allowed limits

**Option 1: Use command-line arguments**
```bash
outerspace pipeline config.toml snakemake_config.yaml \
    --snakemake-args="--executor slurm --jobs 100 \
    --default-resources slurm_account=myaccount slurm_partition=compute mem_mb=4000"
```

**Option 2: Use a profile (recommended)** - Create a profile with `default-resources` as shown in the "Using Profiles" section above. The pipeline will automatically parse and apply these settings to all jobs.

### Workflow Debugging

```bash
# Dry run to check DAG without submitting
outerspace pipeline config.toml snakemake_config.yaml --snakemake-args="--dry-run"

# Print DAG for inspection
outerspace pipeline config.toml snakemake_config.yaml --snakemake-args="--dag | dot -Tpdf > dag.pdf"

# Unlock workflow directory if locked from failed run
outerspace pipeline config.toml snakemake_config.yaml --snakemake-args="--unlock"
```

## Workflow Description

The OUTERSPACE Snakemake workflow (`workflow/Snakefile`) implements a complete analysis pipeline with the following steps:

1. **findseq**: Extracts sequences from FASTQ files based on configuration patterns
   - Supports both single-end and paired-end reads
   - Outputs CSV files with extracted sequences
   - Supports multi-threaded parallel processing (configurable via `threads` parameter)

2. **collapse**: Corrects barcodes using UMI-tools clustering
   - Takes findseq output as input
   - Applies configurable clustering methods
   - Generates corrected barcode files

3. **count**: Counts unique barcodes per key value
   - Processes collapsed barcode files
   - Supports downsampling and filtering
   - Calculates detailed metrics

4. **merge**: Merges multiple csvs into a single file
   - Can do joint UMI correction

5. **stats**: Calculcates sample level metrics about the library
   - See [Stats](docs/stats.md) for more details.


The workflow automatically handles:
- Sample discovery from input directories
- Parallel processing of multiple samples
- Proper file naming and organization
- Dependency tracking between steps

## Workflow Configuration

### Thread Configuration

The `findseq` step supports multi-threaded parallel processing for improved performance on large FASTQ files. You can configure the number of threads in your Snakemake YAML configuration file:

```yaml
# config.yaml
samples: samples.csv
toml: config.toml
threads: 8  # Use 8 threads for findseq processing
```

Alternatively, you can override the thread count when running the pipeline:

```bash
outerspace pipeline config.toml config.yaml --snakemake-args="--config threads=8"
```

Note: The actual speedup depends on your system's I/O performance and CPU capabilities. If not specified, `threads` defaults to 1.

## Snakemake Wrappers

Snakemake wrappers are reusable components that encapsulate tool execution in a standardized way. They are particularly useful for:
- Running analyses on cluster environments
- Ensuring consistent tool execution
- Simplifying workflow development
- Enabling portability across different systems

### OUTERSPACE Wrappers

The OUTERSPACE workflow uses custom wrappers that provide a thin layer around the OUTERSPACE CLI commands. These wrappers:
- Map Snakemake parameters to CLI arguments
- Handle input/output file management
- Provide consistent error handling
- Enable easy integration with cluster schedulers

Each wrapper corresponds to a CLI command:
- `findseq`: Wraps the sequence extraction command
- `collapse`: Wraps the barcode correction command
- `count`: Wraps the barcode counting command
- `merge`: Merges all data into a single csv file.
- `stats`: Calculates sample-level metrics across the inputs.

The wrappers are located in `workflow/wrappers/` and can be customized for specific cluster environments or requirements.

### Using Wrapper Scripts

Wrapper scripts are referenced in the Snakefile using the `script` directive:

```python
rule findseq:
    input:
        reads = get_input_files,
        toml = get_toml_file
    output:
        'findseq/{sample}.csv'
    threads: config.get('threads', 1)
    script:
        'wrappers/findseq/wrapper.py'
```

This approach allows for:
- Easy modification of wrapper behavior
- Consistent execution across environments
- Simple integration with cluster schedulers
- Reproducible analysis workflows


Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.