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
- Handles additional Snakemake arguments
- Provides proper error handling and logging

## Workflow Description

The OUTERSPACE Snakemake workflow (`workflow/Snakefile`) implements a complete analysis pipeline with the following steps:

1. **findseq**: Extracts sequences from FASTQ files based on configuration patterns
   - Supports both single-end and paired-end reads
   - Outputs CSV files with extracted sequences

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

### Using Wrappers

Wrappers are referenced in the Snakefile using the `wrapper` directive:

```python
rule findseq:
    input:
        reads = get_input_files,
        toml = get_toml_file
    output:
        'findseq/{sample}.csv'
    wrapper:
        f'file:{WORKFLOW_DIR}/wrappers/findseq'
```

This approach allows for:
- Easy modification of wrapper behavior
- Consistent execution across environments
- Simple integration with cluster schedulers
- Reproducible analysis workflows


### Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.