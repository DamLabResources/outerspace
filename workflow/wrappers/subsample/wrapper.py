"""Wrapper for outerspace subsample command"""

__author__ = "WND"
__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "0.0.2"

import sys
from outerspace.cli.main import Cli

# This is a common pattern in Snakemake wrappers
if "snakemake" not in locals():
    import snakemake  # type: ignore

# Get input and output files
output_file = snakemake.output[0]
toml_file = snakemake.input.get("toml", None)

# Get CSV input file (exclude the toml file)
if hasattr(snakemake.input, 'csv'):
    # If csv is explicitly named in the input
    csv_file = snakemake.input.csv
else:
    # Filter out the toml file from all inputs
    csv_files = [f for f in snakemake.input if f != toml_file]
    if len(csv_files) != 1:
        raise ValueError(f"Expected exactly one CSV input file, got {len(csv_files)}")
    csv_file = csv_files[0]

# Construct command line arguments
args = [
    'subsample',
    '-c', toml_file,
    '-o', output_file,
]

# Add CSV input file
args.append(csv_file)

# Add threads parameter if available
if hasattr(snakemake, 'threads') and snakemake.threads > 1:
    args.extend(['--threads', str(snakemake.threads)])

# Run the subsample command
cli = Cli(args)
cli.run()

        
# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.

