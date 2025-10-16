"""Barcode counting command for CSV analysis.

This module provides the CountCommand class for counting unique barcodes per key
value in CSV files. It supports both single file and batch processing with
options for downsampling, filtering, and detailed statistics including Gini coefficients.
"""

import csv
import glob
import logging
import os
import random
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Optional, Set, Tuple, Union
from argparse import ArgumentParser

from tqdm import tqdm
import yaml

from outerspace.cli.commands.base import BaseCommand
from outerspace.umi import UMI
from outerspace.stats import GiniCoefficient, SimpsonDiversity
from outerspace.cli.logging_config import setup_logging

# Set up logging
logger = logging.getLogger(__name__)

# Increase CSV field size limit to handle large fields
csv.field_size_limit(sys.maxsize)

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


class CountCommand(BaseCommand):
    """Command for counting unique barcodes per key value in CSV files.

    This command analyzes CSV files to count unique barcodes grouped by key values.
    It provides comprehensive statistics including Gini coefficients for distribution
    analysis and supports various filtering and sampling options.
    """

    def _init_parser(self, subparser: ArgumentParser) -> None:
        """Initialize command-specific argument parser.

        Parameters
        ----------
        subparser : ArgumentParser
            Subparser group to add command arguments to
        """
        parser = subparser.add_parser(
            "count", help="Count unique barcodes per key value in CSV files"
        )

        # Input options (mutually exclusive)
        input_group = parser.add_mutually_exclusive_group(required=True)
        input_group.add_argument("--input-file", help="Input CSV file to process")
        input_group.add_argument(
            "--input-dir", help="Input directory containing CSV files to process"
        )

        # Output options (mutually exclusive)
        output_group = parser.add_mutually_exclusive_group(required=True)
        output_group.add_argument(
            "--output-file", help="Output CSV file for barcode counts"
        )
        output_group.add_argument(
            "--output-dir", help="Output directory for barcode counts"
        )

        # Processing options
        parser.add_argument("--barcode-column", help="Column containing barcodes")
        parser.add_argument("--key-column", help="Column to group by")
        parser.add_argument("--sep", default=",", help="CSV separator (default: ,)")
        parser.add_argument(
            "--row-limit",
            type=int,
            help="Process only the first N rows (for testing)",
            default=None,
        )
        parser.add_argument(
            "--allowed-list",
            help=(
                "DEPRECATED: Use collapse --method allowed/nearest instead. "
            ),
            default=None,
        )
        parser.add_argument(
            "--detailed", action="store_true", help="Include barcode lists in output"
        )
        parser.add_argument(
            "--downsample",
            type=float,
            help="Randomly sample reads with probability between 0 and 1",
        )
        parser.add_argument(
            "--random-seed", type=int, help="Random seed for downsampling"
        )

        self._add_common_args(parser)

    def _read_allowed_keys(self, filepath: str) -> Tuple[str, ...]:
        """Read allowed keys from a text file.

        Parameters
        ----------
        filepath : str
            Path to text file containing allowed keys

        Returns
        -------
        Tuple[str, ...]
            Allowed keys in file order (unique, empties filtered)
        """
        allowed_keys_list = []
        seen: Set[str] = set()
        try:
            with open(filepath, "r") as f:
                for line in f:
                    key = line.strip()
                    if key and key not in seen:
                        seen.add(key)
                        allowed_keys_list.append(key)
            logger.info(f"Loaded {len(allowed_keys_list)} allowed keys from {filepath}")
        except Exception as e:
            logger.error(f"Failed to read allowed keys from {filepath}: {e}")
            raise
        return tuple(allowed_keys_list)

    def _process_single_file(
        self,
        input_file: str,
        output_file: str,
        barcode_col: str,
        key_col: str,
        sep: str,
        row_limit: Optional[int],
        allowed_keys_ordered: Optional[Tuple[str, ...]],
        detailed: bool,
        downsample: Optional[float] = None,
    ) -> Dict[str, Any]:
        """Process a single CSV file and return summary statistics.

        This method reads a CSV file, counts unique barcodes per key, and calculates
        various statistics including Gini coefficients for distribution analysis.

        Parameters
        ----------
        input_file : str
            Path to input CSV file
        output_file : str
            Path to output CSV file
        barcode_col : str
            Column name containing barcodes
        key_col : str
            Column name to group by
        sep : str
            CSV separator character
        row_limit : Optional[int]
            Maximum number of rows to process (for testing)
        allowed_keys_ordered : Optional[Tuple[str, ...]]
            Allowed keys for exact-match filtering only (DEPRECATED - use collapse instead)
        detailed : bool
            Whether to include barcode lists in output
        downsample : Optional[float], default=None
            Probability for random sampling (0-1)

        Returns
        -------
        Dict[str, Any]
            Dictionary containing summary statistics

        Raises
        ------
        ValueError
            If required columns are not found in the input file
        """
        logger.info(f"Processing file: {input_file}")

        # Create UMI objects for statistics calculation
        umi = UMI(mismatches=0)
        key_umi = UMI(mismatches=0)  # For key counts

        allowed_keys_set = set(allowed_keys_ordered) if allowed_keys_ordered else set()

        # Read rows and collect barcodes per key
        barcodes_by_key = defaultdict(set)
        total_rows = 0
        rows_with_allowed_key = 0

        with open(input_file, "r") as f:
            reader = csv.DictReader(f, delimiter=sep)
            headers = reader.fieldnames

            # Verify columns exist
            missing_cols = [col for col in [barcode_col, key_col] if col not in headers]
            if missing_cols:
                raise ValueError(
                    f"Columns not found in input file: {', '.join(missing_cols)}"
                )

            for i, row in enumerate(tqdm(reader, desc="Reading rows")):
                if row_limit and i >= row_limit:
                    break

                # Apply downsampling if specified
                if downsample is not None and random.random() > downsample:
                    continue

                total_rows += 1
                key = str(row[key_col])
                barcode = str(row[barcode_col])

                if not key or not barcode:
                    continue

                # Counting logic
                # Case 1: No allowed_keys list provided -> count all (key, barcode)
                if not allowed_keys_ordered:
                    barcodes_by_key[key].add(barcode)
                    umi.consume(barcode)
                    key_umi.consume(key)
                    continue

                # Case 2: allowed_keys provided and key is in the set (exact match only)
                if key in allowed_keys_set:
                    rows_with_allowed_key += 1
                    barcodes_by_key[key].add(barcode)
                    umi.consume(barcode)
                    key_umi.consume(key)
                    continue

                # Case 3: key not in allowed_keys; ignore (no rescue in count command)
                # Users should use collapse --method nearest for key correction
                continue

        # Calculate summary statistics
        total_keys = len(barcodes_by_key)
        total_barcodes = sum(len(barcodes) for barcodes in barcodes_by_key.values())

        # Log statistics
        logger.info(f"File statistics for {os.path.basename(input_file)}:")
        logger.info(f"Total rows scanned: {total_rows}")
        logger.info(f"Total keys: {total_keys}")
        logger.info(f"Total barcodes: {total_barcodes}")
        logger.info(
            f"Average barcodes per key: {total_barcodes / total_keys if total_keys > 0 else 0:.3f}"
        )

        # Calculate Gini coefficients using the stats module
        barcode_gini = GiniCoefficient.calculate(umi)
        key_gini = GiniCoefficient.calculate(
            key_umi,
            allowed_list=list(allowed_keys_ordered) if allowed_keys_ordered else None,
        )
        logger.info(f"Barcode Gini coefficient: {barcode_gini:.3f}")
        logger.info(f"Key Gini coefficient: {key_gini:.3f}")

        barcode_simpson = SimpsonDiversity.calculate(umi)
        key_simpson = SimpsonDiversity.calculate(
            key_umi,
            allowed_list=list(allowed_keys_ordered) if allowed_keys_ordered else None,
        )
        logger.info(f"Barcode Simpson diversity: {barcode_simpson:.3f}")
        logger.info(f"Key Simpson diversity: {key_simpson:.3f}")

        if allowed_keys_ordered:
            logger.info(f"Rows with allowed key: {rows_with_allowed_key}")
            missing_keys = set(allowed_keys_ordered) - set(barcodes_by_key.keys())
            logger.info(f"Total missing keys: {len(missing_keys)}")
            if len(missing_keys) > 0 and detailed:
                logger.info("Missing keys:")
                for key in sorted(list(missing_keys))[:10]:
                    logger.info(f"  {key}")
                if len(missing_keys) > 10:
                    logger.info(f"  ... and {len(missing_keys) - 10} more")

        # Write output
        self._write_counts(barcodes_by_key, output_file, sep, detailed, key_col)

        return {
            "total_rows": total_rows,
            "total_keys": total_keys,
            "total_barcodes": total_barcodes,
            "barcode_gini": barcode_gini,
            "key_gini": key_gini,
        }

    def _write_counts(
        self,
        barcodes_by_key: Dict[str, Set[str]],
        filepath: str,
        sep: str,
        detailed: bool,
        key_col: str,
    ) -> None:
        """Write barcode counts per key to CSV file.

        Parameters
        ----------
        barcodes_by_key : Dict[str, Set[str]]
            Dictionary mapping keys to sets of barcodes
        filepath : str
            Path to output CSV file
        sep : str
            CSV separator character
        detailed : bool
            Whether to include barcode lists in output
        key_col : str
            Name of the key column
        """
        logger.debug(f"Writing counts to {filepath}")

        with open(filepath, "w", newline="") as f:
            if detailed:
                writer = csv.writer(f, delimiter=sep)
                writer.writerow(
                    [key_col, "unique_barcodes", f"{self.args.barcode_column}_count"]
                )
                for key, barcodes in sorted(barcodes_by_key.items()):
                    writer.writerow([key, ",".join(sorted(barcodes)), len(barcodes)])
            else:
                writer = csv.writer(f, delimiter=sep)
                writer.writerow([key_col, f"{self.args.barcode_column}_count"])
                for key, barcodes in sorted(barcodes_by_key.items()):
                    writer.writerow([key, len(barcodes)])

    def run(self) -> None:
        """Run the count command.

        This method orchestrates the barcode counting process, handling both
        single file and batch processing modes with comprehensive error handling
        and statistics reporting.

        Raises
        ------
        ValueError
            If required arguments are missing or invalid
        """
        # Set up logging
        logger = setup_logging(log_file=self.args.log_file)
        logger.info("Starting barcode counting process")

        # Load config if provided
        if self.args.config:
            self._load_config(self.args.config)

        # Merge config and args with defaults
        defaults = {
            "sep": ",",
            "row_limit": None,
            "detailed": False,
            "downsample": None,
            "random_seed": None,
        }
        self._merge_config_and_args(defaults)

        # Validate required arguments
        if not self.args.barcode_column and not self.args.config:
            raise ValueError("Please provide either --barcode-column or --config")
        if not self.args.key_column and not self.args.config:
            raise ValueError("Please provide either --key-column or --config")

        # Validate input/output arguments
        if not self.args.input_file and not self.args.input_dir:
            raise ValueError("Please provide either --input-file or --input-dir")
        if not self.args.output_file and not self.args.output_dir:
            raise ValueError("Please provide either --output-file or --output-dir")

        # Validate downsampling parameter if provided
        if self.args.downsample is not None:
            if not 0 < self.args.downsample <= 1:
                raise ValueError("Downsample probability must be between 0 and 1")
            if self.args.random_seed is not None:
                random.seed(self.args.random_seed)
                logger.info(f"Using random seed: {self.args.random_seed}")

        # Read allowed keys if specified
        allowed_keys: Optional[Tuple[str, ...]] = None
        if self.args.allowed_list:
            logger.warning("=" * 80)
            logger.warning("DEPRECATION WARNING: --allowed-list in count command")
            logger.warning("")
            logger.warning("The --allowed-list feature is deprecated for the count command.")
            logger.warning("It will only perform exact matching (filtering).")
            logger.warning("")
            logger.warning("For key correction/rescue, use collapse with iterative steps:")
            logger.warning("")
            logger.warning("  [[collapse.steps]]")
            logger.warning("  name = \"key_correction\"")
            logger.warning("  columns = \"your_key_column\"")
            logger.warning("  method = \"nearest\"  # or \"allowed\" for exact matching")
            logger.warning("  allowed_list = \"allowed_keys.txt\"")
            logger.warning("")
            logger.warning("Continuing with exact-match filtering only...")
            logger.warning("=" * 80)
            allowed_keys = self._read_allowed_keys(self.args.allowed_list)

        # No UMI normalization in this command; UMIs are random, we count all

        # Handle single file case
        if self.args.input_file:
            if not os.path.exists(self.args.input_file):
                raise ValueError(f"Input file not found: {self.args.input_file}")

            # Create output directory if needed
            output_path = Path(self.args.output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            try:
                self._process_single_file(
                    input_file=self.args.input_file,
                    output_file=self.args.output_file,
                    barcode_col=self.args.barcode_column,
                    key_col=self.args.key_column,
                    sep=self.args.sep,
                    row_limit=self.args.row_limit,
                    allowed_keys_ordered=allowed_keys,
                    detailed=self.args.detailed,
                    downsample=self.args.downsample,
                )

            except Exception as e:
                logger.error(f"Error processing {self.args.input_file}: {e}")
                raise

            logger.info(
                f"Processing complete. Barcode counts written to: {self.args.output_file}"
            )
            return

        # Handle directory case
        if not os.path.exists(self.args.input_dir):
            raise ValueError(f"Input directory not found: {self.args.input_dir}")

        # Create output directory if it doesn't exist
        os.makedirs(self.args.output_dir, exist_ok=True)

        # Get list of CSV files in input directory
        input_files = glob.glob(os.path.join(self.args.input_dir, "*.csv"))
        if not input_files:
            raise ValueError(f"No CSV files found in {self.args.input_dir}")

        logger.info(f"Found {len(input_files)} CSV files to process")

        # Process each file
        for input_file in tqdm(input_files, desc="Processing files"):
            # Create output filename
            output_file = os.path.join(
                self.args.output_dir, os.path.basename(input_file)
            )

            try:
                self._process_single_file(
                    input_file=input_file,
                    output_file=output_file,
                    barcode_col=self.args.barcode_column,
                    key_col=self.args.key_column,
                    sep=self.args.sep,
                    row_limit=self.args.row_limit,
                    allowed_keys_ordered=allowed_keys,
                    detailed=self.args.detailed,
                    downsample=self.args.downsample,
                )

            except Exception as e:
                logger.error(f"Error processing {input_file}: {e}")
                raise

        logger.info(
            f"Processing complete. Barcode counts written to: {self.args.output_dir}"
        )


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
