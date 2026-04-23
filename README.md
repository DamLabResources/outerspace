# OuterSpace

Outerspace is a collection of tools for analyzing pooled CRISPR screens, viral barcode population studies, and any other application that requires the extraction and couting of variable regions in pooled amplicons.
It contains tools to extract regions of interest, correct sequencing error, assess diversity, and compare between samples.

## Contents
- [Quick Start](#quick-start)
- [Installation](#install)
- [Biologic Applications](docs/usage.md)
- [Basic Commands](docs/commands.md)
- [Configuration](docs/config.md)
- [Detailed Walkthrough](docs/walkthrough.md)
- [Command Line Interface Documentation](docs/cli_api.md)
- [Snakemake Documentation](docs/snakemake.md)

## Quick Start

### Install

```bash
#pip install outerspace
pip install git+https://github.com/DamLabResources/outerspace.git
```

### Design your extraction strategy

`outerspace` uses the `regex` library to extract relevant features from a DNA sequence.
This allows an simple, expressive, and modular strategy for extracting of regions of interest while tolerating mismatches.
It supports both short, paired end reads and long reads.
See the [docs/regex_explainer.md](docs/regex_explainer.md) for a detailed discussion on how to design your extraction strategy.  
[Regex Link](https://pypi.org/project/regex/)

### Create your config file

If you are going to repeating similar experiments often, `outerspace` allows you to encapsulate that information in a `toml` file accepted by all commands.
This ensures repeatability between analyses and can drastically simplify command line execution.
It also facilitates reproducible science as the config can be stored, shared, and tracked.

See the [walkthrough](docs/config.md) for a more detailed discussion on creating your config file.

### Process Your Data

For most analyses, you can use the `pipeline` command to process all your data in one step:

```bash
# Create output directory
mkdir -p results

# Run the pipeline
outerspace pipeline config.toml \
    --input-dir fastq_files \
    --output-dir results \
    --barcode-columns UMI_5prime,UMI_3prime \
    --key-column protospacer \
    --mismatches 2 \
    --method directional \
    --metrics
```

This will:
1. Process all FASTQ files in your input directory
2. Extract sequences using your config patterns
3. Correct barcodes using UMI-tools clustering
4. Count unique barcodes per protospacer
5. Generate metrics files for quality control

For more detailed instructions, including how to run individual commands and perform additional analyses, see the [detailed walkthrough](walkthrough.md).

For running your tasks in parallel or on a cluster consider using our Snakemake [wrappers](wrappers.md).

## Development

The repository expects a **conda environment** in this clone’s `venv` directory (for example
`/home/jupyter-will/repos/outerspace/venv`). Create or refresh it, then use that interpreter
for all tooling:

```bash
# One-time: create env and install with dev + pipeline extras
make venv
# or: conda create -y -p ./venv python=3.12 && ./venv/bin/pip install -e ".[dev,pipeline]"

# Run tools via the env (from repo root)
./venv/bin/python -m tox -e py        # fast unit tests (default)
./venv/bin/python -m tox -e ruff     # lint + format check
./venv/bin/python -m tox -e functional   # long integration/functional tests (optional)

# Or with conda run (same env)
conda run -p /home/jupyter-will/repos/outerspace/venv python -m tox -e py
```

Shorter: `make test` and `make ruff` use `$(CURDIR)/venv` automatically. Install
[pre-commit](https://pre-commit.com/) hooks with `pre-commit install` after the dev extra is installed.

**Documentation (MkDocs):** `pip install -e ".[docs]" && mkdocs serve` (or `tox -e docs`).

- **User guide (Read the Docs):** [https://outerspace.readthedocs.io/](https://outerspace.readthedocs.io/) (enable the project in RTD when ready).

Copyright (C) 2025, SCB, DVK PhD, RB, WND PhD. All rights reserved.
