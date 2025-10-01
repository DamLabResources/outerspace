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
from outerspace.stats import GiniCoefficient
from outerspace.cli.logging_config import setup_logging
from outerspace.nearest import find_closest_umi

# Set up logging
logger = logging.getLogger(__name__)

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
            help="Text file containing allowed keys (one per line)",
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

        # Key rescue options (optional)
        parser.add_argument(
            "--key-rescue",
            action="store_true",
            help=(
                "If --allowed-list is provided, attempt to rescue keys not in the list "
                "by mapping them to the closest allowed key using global alignment."
            ),
        )
        parser.add_argument(
            "--key-min-score",
            type=int,
            default=0,
            help="Minimum alignment score required to rescue a key (default: 0)",
        )
        parser.add_argument(
            "--key-match-score",
            type=int,
            default=1,
            help="Score for matches when aligning keys (default: 1)",
        )
        parser.add_argument(
            "--key-mismatch-penalty",
            type=int,
            default=-1,
            help="Penalty for mismatches when aligning keys (default: -1)",
        )
        parser.add_argument(
            "--key-gap-penalty",
            type=int,
            default=-3,
            help="Penalty for gaps/indels when aligning keys (default: -3)",
        )

        parser.add_argument(
            "--key-rescue-strategy",
            type=str,
            default="random",
            help="Strategy to choose the best mapped key when multiple are found (default: random)",
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
        key_mismatch_penalty: Optional[int],
        key_gap_penalty: Optional[int],
        key_match_score: Optional[int],
        key_min_score: Optional[int],
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
            Allowed keys in file order (for deterministic rescue)
        key_mismatch_penalty : Optional[int]
            Mismatch penalty for key rescue
        key_gap_penalty : Optional[int]
            Gap penalty for key rescue
        key_match_score : Optional[int]
            Match score for key rescue
        key_min_score : Optional[int]
            Minimum score threshold to accept a rescued key
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

        rescued_count = 0
        multiple_rescued_count = 0

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

                # Case 2: allowed_keys provided and key is in the set
                if key in allowed_keys_set:
                    rows_with_allowed_key += 1
                    barcodes_by_key[key].add(barcode)
                    umi.consume(barcode)
                    key_umi.consume(key)
                    continue

                # Case 3: key not in allowed_keys; attempt rescue if enabled
                if self.args.key_rescue:
                    # Map this key to the closest allowed keys; here we use the same
                    # alignment function, passing the set of allowed keys as candidates
                    # Will return a list of keys that are all equally close to the key
                    mapped_keys = find_closest_umi(
                        allowed_keys_ordered,
                        key,
                        mismatch_penalty=key_mismatch_penalty,
                        gap_penalty=key_gap_penalty,
                        match_score=key_match_score,
                        min_score=key_min_score,
                    )
                    if mapped_keys:
                        mapped_key = None
                        # If there is only one mapped key, use it
                        if len(mapped_keys) == 1:
                            mapped_key = mapped_keys[0]
                        else:
                            multiple_rescued_count += 1
                            # If there are multiple mapped keys, use the strategy to choose one
                            if self.args.key_rescue_strategy == "random":
                                mapped_key = random.choice(mapped_keys)
                            elif self.args.key_rescue_strategy == "first":
                                mapped_key = mapped_keys[0]
                            elif self.args.key_rescue_strategy == "last":
                                mapped_key = mapped_keys[-1]
                            elif self.args.key_rescue_strategy == "all":
                                rescued_count += 1
                                for mapped_key in mapped_keys:
                                    barcodes_by_key[mapped_key].add(barcode)
                                    umi.consume(barcode)
                                    key_umi.consume(mapped_key)
                                continue
                            
                        if mapped_key:
                            rescued_count += 1
                            barcodes_by_key[mapped_key].add(barcode)
                            umi.consume(barcode)
                            key_umi.consume(mapped_key)

                        continue
                        
                        
                        

                # If not rescued, ignore this key (do not count)
                continue

        # Calculate summary statistics
        total_keys = len(barcodes_by_key)
        total_barcodes = sum(len(barcodes) for barcodes in barcodes_by_key.values())

        # Calculate Gini coefficients using the stats module
        barcode_result = GiniCoefficient.calculate(umi)
        key_result = GiniCoefficient.calculate(
            key_umi,
            allowed_list=list(allowed_keys_ordered) if allowed_keys_ordered else None,
        )
        barcode_gini = barcode_result
        key_gini = key_result

        # Log statistics
        logger.info(f"File statistics for {os.path.basename(input_file)}:")
        logger.info(f"Total rows scanned: {total_rows}")
        logger.info(f"Total keys: {total_keys}")
        logger.info(f"Total barcodes: {total_barcodes}")
        logger.info(
            f"Average barcodes per key: {total_barcodes / total_keys if total_keys > 0 else 0:.3f}"
        )
        logger.info(f"Barcode Gini coefficient: {barcode_gini:.3f}")
        logger.info(f"Key Gini coefficient: {key_gini:.3f}")

        # Log rescue info if applicable
        if self.args.key_rescue and allowed_keys_ordered:
            logger.info(f"Keys rescued (mapped to allowed list): {rescued_count}")
            logger.info(f"Keys rescued (mapped to allowed list with multiple matches): {multiple_rescued_count}")

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
                    key_mismatch_penalty=self.args.key_mismatch_penalty,
                    key_gap_penalty=self.args.key_gap_penalty,
                    key_match_score=self.args.key_match_score,
                    key_min_score=self.args.key_min_score,
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
                    key_mismatch_penalty=self.args.key_mismatch_penalty,
                    key_gap_penalty=self.args.key_gap_penalty,
                    key_match_score=self.args.key_match_score,
                    key_min_score=self.args.key_min_score,
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
