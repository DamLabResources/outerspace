"""Wrapper for outerspace collapse barcode correction"""

__author__ = "WND"
__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__email__ = "wnd22@drexel.edu"
__license__ = "MIT"
__version__ = "0.0.2"

from outerspace.cli.main import Cli

# This is a common pattern in Snakemake wrappers
if "snakemake" not in locals():
    import snakemake  # type: ignore

# Get input and output files
input_file = snakemake.input[0]
toml_file = snakemake.input.get("toml", None)
output_file = snakemake.output[0]


# Construct command line arguments
args = [
    'collapse',
    '--input-file', input_file,
    '--output-file', output_file,
    '-c', toml_file,
    '--threads', str(snakemake.threads),
]

# Run the collapse command
cli = Cli(args)
cli.run() 

# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.