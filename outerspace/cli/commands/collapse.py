"""Barcode correction command for UMI clustering.

This module provides the CollapseCommand class for correcting barcodes in CSV files
using UMI-tools clustering. It supports both single file and batch processing
with various clustering methods and comprehensive metrics reporting.
"""

import csv
import glob
import logging
import os
import random
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple
from argparse import ArgumentParser

from tqdm import tqdm

from outerspace.cli.commands.base import BaseCommand
from outerspace.umi import UMI
from outerspace.nearest import NearestUMIFinder

# Set up logging
logger = logging.getLogger(__name__)

# Increase CSV field size limit to handle large fields
csv.field_size_limit(sys.maxsize)

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


class CollapseCommand(BaseCommand):
    """Command for correcting barcodes in CSV files using UMI clustering.

    This command performs barcode correction using UMI-tools clustering algorithms.
    It supports processing single files or entire directories, with options for
    different clustering methods and comprehensive metrics reporting.
    """

    def _init_parser(self, subparser: ArgumentParser) -> None:
        """Initialize command-specific argument parser.

        Parameters
        ----------
        subparser : ArgumentParser
            Subparser group to add command arguments to
        """
        parser = subparser.add_parser(
            "collapse", help="Correct barcodes in CSV files using UMI-tools clustering"
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
            "--output-file", help="Output CSV file for corrected barcodes"
        )
        output_group.add_argument(
            "--output-dir", help="Output directory for corrected CSV files"
        )

        # Processing options
        parser.add_argument(
            "--columns",
            help="Column(s) containing barcodes to correct. Can be a single column or comma-separated list",
        )
        parser.add_argument(
            "--mismatches",
            type=int,
            default=2,
            help="Number of mismatches allowed for clustering (default: 2)",
        )
        parser.add_argument("--sep", default=",", help="CSV separator (default: ,)")
        parser.add_argument(
            "--row-limit",
            type=int,
            help="Process only the first N rows (for testing)",
            default=None,
        )
        parser.add_argument(
            "--method",
            choices=["cluster", "adjacency", "directional", "allowed", "nearest"],
            default="directional",
            help=(
                "Correction method: 'cluster', 'adjacency', 'directional' for UMI clustering; "
                "'allowed' for exact matching with --allowed-list; "
                "'nearest' for nearest-neighbor rescue with --allowed-list (default: directional)"
            ),
        )
        parser.add_argument("--metrics", help="Output YAML file for metrics")
        
        # Allowed list options (required for method='allowed' or 'nearest')
        parser.add_argument(
            "--allowed-list",
            help=(
                "Text file containing allowed values for --columns (one per line). "
                "Required for --method allowed or --method nearest."
            ),
            default=None,
        )
        # Alignment scoring parameters (for use with --method nearest)
        parser.add_argument(
            "--min-score",
            type=int,
            default=0,
            help="Minimum alignment score required to rescue a value (default: 0). For use with --method nearest.",
        )
        parser.add_argument(
            "--match-score",
            type=int,
            default=1,
            help="Score for matches when aligning (default: 1). For use with --method nearest.",
        )
        parser.add_argument(
            "--mismatch-penalty",
            type=int,
            default=-1,
            help="Penalty for mismatches when aligning (default: -1). For use with --method nearest.",
        )
        parser.add_argument(
            "--gap-penalty",
            type=int,
            default=-3,
            help="Penalty for gaps/indels when aligning (default: -3). For use with --method nearest.",
        )
        
        # K-mer prescreen parameters (for use with --method nearest)
        parser.add_argument(
            "--rescue-kmer-size",
            type=int,
            default=3,
            help="K-mer size for approximate prescreen (default: 3). For use with --method nearest.",
        )
        parser.add_argument(
            "--rescue-min-overlap",
            type=int,
            default=1,
            help=(
                "Minimum number of shared k-mers required before alignment "
                "(default: 1). For use with --method nearest."
            ),
        )
        parser.add_argument(
            "--rescue-exhaustive",
            action="store_true",
            default=False,
            help=(
                "Disable k-mer prescreen and align against all allowed values "
                "(slower, exhaustive search). For use with --method nearest."
            ),
        )
        parser.add_argument(
            "--rescue-strategy",
            type=str,
            default="random",
            help="Strategy to choose the best match when multiple are found (default: random). For use with --method nearest.",
        )
        
        # Threading options
        parser.add_argument(
            "--threads",
            type=int,
            default=1,
            help="Number of threads for parallel processing (default: 1). For use with --method nearest.",
        )
        
        self._add_common_args(parser)

    def _parse_columns(self, columns_str: str) -> List[str]:
        """Parse comma-separated column string into list of column names.

        Parameters
        ----------
        columns_str : str
            Comma-separated string of column names

        Returns
        -------
        List[str]
            List of column names with whitespace stripped
        """
        return [col.strip() for col in columns_str.split(",")]

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
        columns: List[str],
        mismatches: int,
        sep: str,
        row_limit: Optional[int],
        method: str,
        allowed_keys_ordered: Optional[Tuple[str, ...]] = None,
        mismatch_penalty: Optional[int] = None,
        gap_penalty: Optional[int] = None,
        match_score: Optional[int] = None,
        min_score: Optional[int] = None,
        rescue_exhaustive: Optional[bool] = None,
        rescue_kmer_size: Optional[int] = None,
        rescue_min_overlap: Optional[int] = None,
        rescue_strategy: Optional[str] = None,
        threads: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Process a single CSV file and return metrics.

        This method reads a CSV file, performs correction using either UMI clustering
        (for methods: cluster, adjacency, directional), exact matching (for method: allowed),
        or nearest-neighbor matching (for method: nearest). It writes the corrected data
        to an output file.

        Parameters
        ----------
        input_file : str
            Path to input CSV file
        output_file : str
            Path to output CSV file
        columns : List[str]
            List of column names to correct
        mismatches : int
            Number of mismatches allowed for UMI clustering (not used for method='allowed' or 'nearest')
        sep : str
            CSV separator character
        row_limit : Optional[int]
            Maximum number of rows to process (for testing)
        method : str
            Correction method: 'cluster', 'adjacency', 'directional' for UMI clustering,
            'allowed' for exact matching with allowed list, or 'nearest' for nearest-neighbor
            matching with allowed list
        allowed_keys_ordered : Optional[Tuple[str, ...]], default=None
            Allowed values in file order (required for method='allowed' or 'nearest')
        mismatch_penalty : Optional[int], default=None
            Mismatch penalty for nearest-neighbor alignment
        gap_penalty : Optional[int], default=None
            Gap penalty for nearest-neighbor alignment
        match_score : Optional[int], default=None
            Match score for nearest-neighbor alignment
        min_score : Optional[int], default=None
            Minimum score threshold to accept a match
        rescue_exhaustive : Optional[bool], default=None
            Whether to disable k-mer prescreen (exhaustive search)
        rescue_kmer_size : Optional[int], default=None
            K-mer size for prescreening
        rescue_min_overlap : Optional[int], default=None
            Minimum k-mer overlap for prescreening
        rescue_strategy : Optional[str], default=None
            Strategy for choosing among multiple matches
        threads : Optional[int], default=None
            Number of threads for parallel processing (for method='nearest')

        Returns
        -------
        Dict[str, Any]
            Dictionary containing processing metrics

        Raises
        ------
        ValueError
            If required columns are not found in the input file
        """
        logger.info(f"Processing file: {input_file}")

        # Read all rows first
        rows = []
        with open(input_file, "r") as f:
            reader = csv.DictReader(f, delimiter=sep)
            headers = reader.fieldnames

            # Verify columns exist
            missing_cols = [col for col in columns if col not in headers]
            if missing_cols:
                raise ValueError(
                    f"Columns not found in input file: {', '.join(missing_cols)}"
                )

            for i, row in enumerate(tqdm(reader, desc="Reading rows")):
                if row_limit and i >= row_limit:
                    break
                rows.append(row)

        if row_limit:
            logger.info(f"Processing first {row_limit} rows of {input_file}")

        # Determine corrected column name
        corrected_col_name = "_".join(columns) + "_corrected"

        # Process based on method
        if method in ["allowed", "nearest"]:
            # Use allowed list matching (exact for 'allowed', nearest-neighbor for 'nearest')
            if allowed_keys_ordered is None or len(allowed_keys_ordered) == 0:
                raise ValueError(f"Method '{method}' requires a non-empty --allowed-list")
            
            if len(columns) > 1:
                raise ValueError(f"Method '{method}' currently only supports single column correction")
            
            column = columns[0]
            allowed_keys_set = set(allowed_keys_ordered)
            rescued_count = 0
            multiple_rescued_count = 0
            
            # Prepare nearest finder for 'nearest' method
            finder = None
            if method == "nearest":
                finder = NearestUMIFinder(
                    allowed_list=allowed_keys_ordered,
                    mismatch_penalty=mismatch_penalty if mismatch_penalty is not None else -1,
                    gap_penalty=gap_penalty if gap_penalty is not None else -3,
                    match_score=match_score if match_score is not None else 1,
                    min_score=min_score if min_score is not None else 0,
                    use_prescreen=not (rescue_exhaustive if rescue_exhaustive is not None else False),
                    kmer_size=rescue_kmer_size if rescue_kmer_size is not None else 3,
                    min_kmer_overlap=rescue_min_overlap if rescue_min_overlap is not None else 1,
                )
            
            corrected_rows = []
            
            if method == "nearest" and finder:
                # Use batch processing for nearest-neighbor matching
                
                # First pass: identify exact matches and collect rescue candidates
                logger.info("Scanning rows for rescue candidates...")
                rescue_indices = []
                rescue_values = []
                
                for i, row in enumerate(rows):
                    value = str(row.get(column, ""))
                    if not value or value in allowed_keys_set:
                        # Exact match or empty - no rescue needed
                        pass
                    else:
                        # Need rescue
                        rescue_indices.append(i)
                        rescue_values.append(value)
                
                logger.info(f"Found {len(rescue_values)} unique values requiring nearest-neighbor matching")
                
                # Batch process rescue values using find_many (with progress bar inside)
                rescue_results = {}
                if rescue_values:
                    num_threads = threads if threads is not None else 1
                    logger.info(f"Using {num_threads} thread(s) for parallel processing")
                    # Get unique rescue values to avoid duplicate work
                    unique_rescue_values = list(set(rescue_values))
                    logger.info(f"Processing {len(unique_rescue_values)} unique rescue candidates")
                    results = finder.find_many(unique_rescue_values, threads=num_threads)
                    rescue_results = dict(zip(unique_rescue_values, results))
                
                # Second pass: apply corrections
                logger.info("Applying corrections to all rows...")
                for row in rows:
                    corrected_row = row.copy()
                    value = str(row.get(column, ""))
                    
                    if not value:
                        corrected_row[corrected_col_name] = ""
                    elif value in allowed_keys_set:
                        # Exact match
                        corrected_row[corrected_col_name] = value
                    else:
                        # Check rescue results
                        mapped_values = rescue_results.get(value)
                        if mapped_values:
                            if len(mapped_values) == 1:
                                corrected_row[corrected_col_name] = mapped_values[0]
                                rescued_count += 1
                            else:
                                multiple_rescued_count += 1
                                # Use strategy to choose one
                                strategy = rescue_strategy if rescue_strategy is not None else "random"
                                if strategy == "random":
                                    corrected_row[corrected_col_name] = random.choice(mapped_values)
                                elif strategy == "first":
                                    corrected_row[corrected_col_name] = mapped_values[0]
                                elif strategy == "last":
                                    corrected_row[corrected_col_name] = mapped_values[-1]
                                else:
                                    corrected_row[corrected_col_name] = mapped_values[0]
                                rescued_count += 1
                        else:
                            # No rescue found, set to empty
                            corrected_row[corrected_col_name] = ""
                    
                    corrected_rows.append(corrected_row)
            else:
                # Method is 'allowed' - use simple exact matching
                desc = "Exact matching with allowed list"
                for row in tqdm(rows, desc=desc):
                    corrected_row = row.copy()
                    value = str(row.get(column, ""))
                    
                    if not value:
                        corrected_row[corrected_col_name] = ""
                    elif value in allowed_keys_set:
                        # Exact match
                        corrected_row[corrected_col_name] = value
                    else:
                        # method=='allowed' and not in allowed list, set to empty
                        corrected_row[corrected_col_name] = ""
                    
                    corrected_rows.append(corrected_row)
            
            if method == "nearest":
                logger.info(f"Values rescued (mapped to allowed list): {rescued_count}")
                logger.info(f"Values rescued with multiple matches: {multiple_rescued_count}")
            
            # Generate metrics
            unique_before = len(set(str(row.get(column, "")) for row in rows if row.get(column, "")))
            unique_after = len(set(row[corrected_col_name] for row in corrected_rows if row[corrected_col_name]))
            
            metrics = {
                "correction_stats": {
                    "unique_values_before": unique_before,
                    "unique_values_after": unique_after,
                    "values_rescued": rescued_count if method == "nearest" else 0,
                },
            }
            
        else:
            # Use UMI clustering (cluster, adjacency, directional)
            umi = UMI(mismatches=mismatches, method=method)

            # Add barcodes to UMI object
            for row in rows:
                # Join multiple columns if specified
                if len(columns) > 1:
                    # Check if all columns are present
                    if all(row.get(col, "") for col in columns):
                        combined_bc = "".join(str(row[col]) for col in columns)
                        umi.consume(combined_bc)
                    else:
                        logger.debug(f"Skipping row with missing columns: {row}")
                else:
                    combined_bc = str(row[columns[0]])
                    if combined_bc:  # Skip empty values
                        umi.consume(combined_bc)

            # Create mapping
            logger.info(
                f"Creating clusters from {len(umi._counts)} unique barcodes from {len(rows)} rows "
                f"with {mismatches} mismatches using {method} method"
            )
            umi.create_mapping()

            # Correct barcodes in rows
            corrected_rows = []
            for row in tqdm(rows, desc="Correcting barcodes"):
                corrected_row = row.copy()

                # Join multiple columns if specified
                if len(columns) > 1:
                    combined_bc = "".join(str(row[col]) for col in columns)
                else:
                    combined_bc = str(row[columns[0]])

                # Create corrected column
                if combined_bc:
                    try:
                        corrected = umi[combined_bc]
                        corrected_row[corrected_col_name] = corrected.decode("ascii")
                    except KeyError:
                        corrected_row[corrected_col_name] = combined_bc
                else:
                    corrected_row[corrected_col_name] = ""

                corrected_rows.append(corrected_row)

            # Generate metrics for UMI clustering
            metrics = {
                "barcode_counts": {
                    "unique_barcodes_before": len(umi._counts),
                    "unique_barcodes_after": len(umi.corrected_counts),
                    "total_reads": sum(umi._counts.values()),
                },
                "correction_details": {
                    "clusters_formed": len(set(umi._mapping.values())),
                    "barcodes_corrected": len(umi._mapping) - len(umi.corrected_counts),
                },
            }

        # Write output
        self._write_csv(corrected_rows, output_file, sep)

        logger.info(f"Processed {len(rows)} rows using {method} method")
        return metrics

    def _write_csv(self, rows: List[Dict[str, str]], filepath: str, sep: str) -> None:
        """Write rows to CSV file.

        Parameters
        ----------
        rows : List[Dict[str, str]]
            List of dictionaries representing CSV rows
        filepath : str
            Path to output CSV file
        sep : str
            CSV separator character
        """
        if not rows:
            logger.warning("No rows to write")
            return

        logger.debug(f"Writing {len(rows)} rows to {filepath}")

        with open(filepath, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter=sep)
            writer.writeheader()
            writer.writerows(rows)

    def _write_metrics(self, metrics: Dict[str, Any], filepath: str) -> None:
        """Write metrics to YAML file.

        Parameters
        ----------
        metrics : Dict[str, Any]
            Dictionary of metrics to write
        filepath : str
            Path to output YAML file
        """
        import yaml

        logger.info(f"Writing metrics to {filepath}")

        with open(filepath, "w") as f:
            yaml.dump(metrics, f, default_flow_style=False)

    def _execute_steps(
        self,
        steps: List[Dict[str, Any]],
        input_path: str,
        output_path: str,
        is_single_file: bool = False,
    ) -> None:
        """Execute a series of collapse steps iteratively.
        
        Parameters
        ----------
        steps : List[Dict[str, Any]]
            List of step configurations from the config file
        input_path : str
            Initial input file or directory
        output_path : str
            Final output file or directory
        is_single_file : bool, optional
            If True, process a single file; if False, process a directory
        """
        mode = "single file" if is_single_file else "batch"
        logger.info(f"Executing {len(steps)} collapse steps iteratively ({mode} mode)")
        
        temp_paths = []
        current_input = input_path
        
        try:
            for i, step in enumerate(steps):
                step_num = i + 1
                is_last_step = (step_num == len(steps))
                
                logger.info(f"=== Step {step_num}/{len(steps)}: {step.get('name', f'step_{step_num}')} ===")
                
                # Determine output path
                if is_last_step:
                    current_output = output_path
                else:
                    # Create temp path for intermediate results
                    if is_single_file:
                        # For single file, create a temp file
                        fd, temp_file = tempfile.mkstemp(suffix=".csv", prefix=f"collapse_step{step_num}_")
                        os.close(fd)  # Close the file descriptor
                        temp_paths.append(temp_file)
                        current_output = temp_file
                    else:
                        # For directory, create a temp directory
                        temp_dir = tempfile.mkdtemp(prefix=f"collapse_step{step_num}_")
                        temp_paths.append(temp_dir)
                        current_output = temp_dir
                
                # Extract step parameters
                columns_str = step.get("columns")
                if not columns_str:
                    raise ValueError(f"Step {step_num}: 'columns' is required")
                columns = self._parse_columns(columns_str)
                
                method = step.get("method", "directional")
                mismatches = step.get("mismatches", 2)
                sep = step.get("sep", ",")
                
                # Read allowed list if specified
                allowed_keys: Optional[Tuple[str, ...]] = None
                if "allowed_list" in step:
                    allowed_list_path = step["allowed_list"]
                    # Handle relative paths relative to config file
                    if self.args.config and not os.path.isabs(allowed_list_path):
                        config_dir = os.path.dirname(self.args.config)
                        allowed_list_path = os.path.join(config_dir, allowed_list_path)
                    allowed_keys = self._read_allowed_keys(allowed_list_path)
                
                # Validate method requirements
                if method in ["allowed", "nearest"]:
                    if not allowed_keys or len(allowed_keys) == 0:
                        raise ValueError(f"Step {step_num}: Method '{method}' requires 'allowed_list'")
                
                if is_single_file:
                    # Single file mode
                    logger.info(f"Processing file: {os.path.basename(current_input)} → {os.path.basename(current_output)}")
                    
                    # Create output directory if needed
                    if is_last_step:
                        output_dir = os.path.dirname(current_output)
                        if output_dir:
                            os.makedirs(output_dir, exist_ok=True)
                    
                    self._process_single_file(
                        current_input,
                        current_output,
                        columns,
                        mismatches,
                        sep,
                        None,  # row_limit
                        method,
                        allowed_keys_ordered=allowed_keys,
                        mismatch_penalty=step.get("mismatch_penalty", -1),
                        gap_penalty=step.get("gap_penalty", -3),
                        match_score=step.get("match_score", 1),
                        min_score=step.get("min_score", 0),
                        rescue_exhaustive=step.get("rescue_exhaustive", False),
                        rescue_kmer_size=step.get("rescue_kmer_size", 3),
                        rescue_min_overlap=step.get("rescue_min_overlap", 1),
                        rescue_strategy=step.get("rescue_strategy", "random"),
                        threads=step.get("threads", self.args.threads if hasattr(self.args, "threads") else 1),
                    )
                else:
                    # Directory mode
                    # Create output directory
                    os.makedirs(current_output, exist_ok=True)
                    
                    # Get list of files to process
                    input_files = glob.glob(os.path.join(current_input, "*.csv"))
                    if not input_files:
                        raise ValueError(f"Step {step_num}: No CSV files found in {current_input}")
                    
                    logger.info(f"Processing {len(input_files)} files: {current_input} → {current_output}")
                    
                    # Process each file
                    for input_file in tqdm(input_files, desc=f"Step {step_num}"):
                        output_file = os.path.join(current_output, os.path.basename(input_file))
                        
                        self._process_single_file(
                            input_file,
                            output_file,
                            columns,
                            mismatches,
                            sep,
                            None,  # row_limit
                            method,
                            allowed_keys_ordered=allowed_keys,
                            mismatch_penalty=step.get("mismatch_penalty", -1),
                            gap_penalty=step.get("gap_penalty", -3),
                            match_score=step.get("match_score", 1),
                            min_score=step.get("min_score", 0),
                            rescue_exhaustive=step.get("rescue_exhaustive", False),
                            rescue_kmer_size=step.get("rescue_kmer_size", 3),
                            rescue_min_overlap=step.get("rescue_min_overlap", 1),
                            rescue_strategy=step.get("rescue_strategy", "random"),
                            threads=step.get("threads", self.args.threads if hasattr(self.args, "threads") else 1),
                        )
                
                logger.info(f"Step {step_num} complete")
                
                # Update input for next iteration
                current_input = current_output
        
        finally:
            # Clean up temporary paths
            for temp_path in temp_paths:
                if os.path.exists(temp_path):
                    if os.path.isfile(temp_path):
                        os.remove(temp_path)
                        logger.debug(f"Cleaned up temp file: {temp_path}")
                    elif os.path.isdir(temp_path):
                        shutil.rmtree(temp_path)
                        logger.debug(f"Cleaned up temp directory: {temp_path}")

    def run(self) -> None:
        """Run the collapse command.

        This method orchestrates the barcode correction process, handling both
        single file and batch processing modes with comprehensive error handling
        and metrics reporting.

        Raises
        ------
        ValueError
            If required arguments are missing or invalid
        """
        logger.info("Starting barcode correction process")

        # Load config if provided
        if self.args.config:
            self._load_config(self.args.config)

        # Check if steps are defined in config
        if self._config and "steps" in self._config and isinstance(self._config["steps"], list):
            # Execute steps mode
            logger.info("Running in steps mode (iterative collapse)")
            
            # Validate input/output arguments
            if not (self.args.input_file or self.args.input_dir):
                raise ValueError("Steps mode requires either --input-file or --input-dir")
            if not (self.args.output_file or self.args.output_dir):
                raise ValueError("Steps mode requires either --output-file or --output-dir")
            
            # Determine if we're in single file or directory mode
            if self.args.input_file:
                if not self.args.output_file:
                    raise ValueError("Steps mode with --input-file requires --output-file")
                
                self._execute_steps(
                    self._config["steps"],
                    self.args.input_file,
                    self.args.output_file,
                    is_single_file=True,
                )
                logger.info(f"All steps complete. Final output: {self.args.output_file}")
            else:
                if not self.args.output_dir:
                    raise ValueError("Steps mode with --input-dir requires --output-dir")
                
                self._execute_steps(
                    self._config["steps"],
                    self.args.input_dir,
                    self.args.output_dir,
                    is_single_file=False,
                )
                logger.info(f"All steps complete. Final output in: {self.args.output_dir}")
            
            return

        # Merge config and args with defaults
        defaults = {
            "mismatches": 2,
            "sep": ",",
            "method": "directional",
            "row_limit": None,
            "threads": 1,
        }
        self._merge_config_and_args(defaults)

        # Validate required arguments
        if not self.args.columns and not self.args.config:
            raise ValueError("Please provide either --columns or --config")

        # Parse columns argument
        columns = self._parse_columns(self.args.columns) if self.args.columns else []

        # Read allowed keys if specified
        allowed_keys: Optional[Tuple[str, ...]] = None
        if self.args.allowed_list:
            allowed_keys = self._read_allowed_keys(self.args.allowed_list)
        
        # Validate parameters for method='allowed' or 'nearest'
        if self.args.method in ["allowed", "nearest"] and not self.args.allowed_list:
            raise ValueError(f"Method '{self.args.method}' requires --allowed-list to be specified")

        # Validate input/output arguments
        if not self.args.input_file and not self.args.input_dir:
            raise ValueError("Please provide either --input-file or --input-dir")
        if not self.args.output_file and not self.args.output_dir:
            raise ValueError("Please provide either --output-file or --output-dir")

        # Handle single file case
        if self.args.input_file:
            if not self.args.output_file:
                raise ValueError(
                    "Please provide an output file when using --input-file"
                )

            # Validate input file exists
            if not os.path.exists(self.args.input_file):
                raise ValueError(f"Input file not found: {self.args.input_file}")

            # Create output directory if needed
            output_path = Path(self.args.output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            try:
                metrics = self._process_single_file(
                    self.args.input_file,
                    self.args.output_file,
                    columns,
                    self.args.mismatches,
                    self.args.sep,
                    self.args.row_limit,
                    self.args.method,
                    allowed_keys_ordered=allowed_keys,
                    mismatch_penalty=self.args.mismatch_penalty,
                    gap_penalty=self.args.gap_penalty,
                    match_score=self.args.match_score,
                    min_score=self.args.min_score,
                    rescue_exhaustive=self.args.rescue_exhaustive,
                    rescue_kmer_size=self.args.rescue_kmer_size,
                    rescue_min_overlap=self.args.rescue_min_overlap,
                    rescue_strategy=self.args.rescue_strategy,
                    threads=self.args.threads,
                )

                # Print metrics
                logger.info(f"Metrics for {os.path.basename(self.args.input_file)}:")
                for category, values in metrics.items():
                    logger.info(f"{category}:")
                    for key, value in values.items():
                        logger.info(f"  {key}: {value}")

                # Write metrics to YAML file if specified
                if self.args.metrics:
                    self._write_metrics(
                        {os.path.basename(self.args.input_file): metrics},
                        self.args.metrics,
                    )
                    logger.info(f"Metrics written to: {self.args.metrics}")

            except Exception as e:
                logger.error(f"Error processing {self.args.input_file}: {e}")
                raise

            logger.info(
                f"Processing complete. Corrected file written to: {self.args.output_file}"
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

        # Collect metrics for all files
        all_metrics = {}

        # Process each file
        for input_file in tqdm(input_files, desc="Processing files"):
            # Create output filename
            output_file = os.path.join(
                self.args.output_dir, os.path.basename(input_file)
            )

            try:
                metrics = self._process_single_file(
                    input_file,
                    output_file,
                    columns,
                    self.args.mismatches,
                    self.args.sep,
                    self.args.row_limit,
                    self.args.method,
                    allowed_keys_ordered=allowed_keys,
                    mismatch_penalty=self.args.mismatch_penalty,
                    gap_penalty=self.args.gap_penalty,
                    match_score=self.args.match_score,
                    min_score=self.args.min_score,
                    rescue_exhaustive=self.args.rescue_exhaustive,
                    rescue_kmer_size=self.args.rescue_kmer_size,
                    rescue_min_overlap=self.args.rescue_min_overlap,
                    rescue_strategy=self.args.rescue_strategy,
                    threads=self.args.threads,
                )

                # Store metrics for this file
                all_metrics[os.path.basename(input_file)] = metrics

                # Print metrics for this file
                logger.info(f"Metrics for {os.path.basename(input_file)}:")
                for category, values in metrics.items():
                    logger.info(f"{category}:")
                    for key, value in values.items():
                        logger.info(f"  {key}: {value}")

            except Exception as e:
                logger.error(f"Error processing {input_file}: {e}")
                raise

        # Write metrics to YAML file if specified
        if self.args.metrics:
            self._write_metrics(all_metrics, self.args.metrics)
            logger.info(f"Metrics written to: {self.args.metrics}")

        logger.info(
            f"Processing complete. Corrected files written to: {self.args.output_dir}"
        )


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
