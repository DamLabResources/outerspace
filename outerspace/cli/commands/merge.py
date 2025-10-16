"""Merge command for combining multiple UMI count files.

This module provides the MergeCommand class for merging multiple UMI count files
into a single file. It supports various output formats.
"""

import logging
import os
from pathlib import Path
from argparse import ArgumentParser

from outerspace.cli.commands.base import BaseCommand
from outerspace.umi import UmiCollection
from outerspace.cli.logging_config import setup_logging

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


class MergeCommand(BaseCommand):
    """Command for merging multiple UMI count files.

    This command combines multiple CSV files containing UMI counts into a single
    output file. It supports both wide and long output formats.
    """

    def _init_parser(self, subparser: ArgumentParser) -> None:
        """Initialize command-specific argument parser.

        Parameters
        ----------
        subparser : ArgumentParser
            Subparser group to add command arguments to
        """
        parser = subparser.add_parser(
            "merge", help="Merge multiple UMI count files into a single file"
        )
        parser.add_argument("files", nargs="+", help="Input CSV files to merge")
        parser.add_argument(
            "--output-file", required=True, help="Output CSV file for merged counts"
        )
        parser.add_argument("--key-column", help="Column containing UMIs")
        parser.add_argument(
            "--count-column",
            help="Column containing counts (if not provided, assumes count=1)",
        )
        parser.add_argument(
            "--sample-names",
            nargs="+",
            help="Optional list of sample names (must match number of input files)",
        )
        parser.add_argument("--sep", default=",", help="CSV separator (default: ,)")
        parser.add_argument(
            "--format",
            choices=["wide", "long"],
            default="wide",
            help="Output format: wide (samples as columns) or long (sample,umi,count columns)",
        )
        self._add_common_args(parser)

    def run(self) -> None:
        """Run the merge command.

        This method orchestrates the file merging process, handling configuration
        loading, UMI collection creation, and output generation with comprehensive
        error handling.

        Raises
        ------
        ValueError
            If required arguments are missing or invalid
        """
        # Set up logging
        logger = setup_logging(log_file=self.args.log_file)
        logger.info("Starting UMI file merge process")

        # Load config if provided
        if self.args.config:
            self._load_config(self.args.config)

        # Merge config and args with defaults
        defaults = {
            "sep": ",",
            "sample_names": None,
            "format": "wide",
        }
        self._merge_config_and_args(defaults)

        # Validate required arguments
        if not self.args.key_column and not self.args.config:
            raise ValueError("Please provide either --key-column or --config")

        # Validate input files exist
        for file in self.args.files:
            if not os.path.exists(file):
                raise ValueError(f"Input file not found: {file}")

        # Create output directory if needed
        output_path = Path(self.args.output_file)
        if output_path.parent:  # Only create directory if there is a path component
            output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            logger.info(f"Starting merge of {len(self.args.files)} files")

            # Create UmiCollection from input files
            logger.info("Creating UmiCollection from input files")
            collection = UmiCollection.from_csvs(
                self.args.files,
                column=self.args.key_column,
                sample_names=self.args.sample_names,
                sep=self.args.sep,
                count_column=self.args.count_column,
            )

            # Write merged data to output file
            logger.info(f"Writing merged data to {self.args.output_file}")
            collection.write(
                self.args.output_file, 
                sep=self.args.sep, 
                format=self.args.format,
                key_name=self.args.key_column
            )

            logger.info(
                f"Successfully merged {len(self.args.files)} files into: {self.args.output_file}"
            )

        except Exception as e:
            logger.error(f"Error merging files: {e}")
            raise


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
