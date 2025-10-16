"""Subsample command for estimating metric stability across sample sizes.

This module provides the SubsampleCommand class for performing random subsampling
on collapse files at various sample sizes to estimate metric stability.
"""

import csv
import logging
import os
import sys
import tempfile
from argparse import ArgumentParser
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
from tqdm import tqdm

from outerspace.cli.commands.base import BaseCommand
from outerspace.stats import (
    BaseStatistic,
    GiniCoefficient,
    ShannonDiversity,
    SimpsonDiversity,
    HillNumber,
    UMIRecoveryRate,
    UMIEfficiencyRate,
    ErrorRate,
)

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

# Registry mapping metric method names to statistic classes
STAT_REGISTRY: Dict[str, type[BaseStatistic]] = {
    "gini_coefficient": GiniCoefficient,
    "shannon_diversity": ShannonDiversity,
    "simpson_diversity": SimpsonDiversity,
    "hill_number": HillNumber,
    "umi_recovery_rate": UMIRecoveryRate,
    "umi_efficiency_rate": UMIEfficiencyRate,
    "error_rate": ErrorRate,
}


class SubsampleCommand(BaseCommand):
    """Command for estimating metric stability through subsampling.

    This command performs random subsampling on a collapse file at various
    sample sizes with multiple replicates to estimate metric stability.
    """

    def _init_parser(self, subparser: ArgumentParser) -> None:
        """Initialize command-specific argument parser.

        Parameters
        ----------
        subparser : ArgumentParser
            Subparser group to add command arguments to
        """
        parser = subparser.add_parser(
            "subsample",
            help="Estimate metric stability through random subsampling",
        )
        parser.add_argument(
            "input_file",
            help="Input CSV file to subsample (collapse output)",
        )
        parser.add_argument("--sep", default=",", help="CSV separator (default: ,)")
        parser.add_argument(
            "--sample-sizes",
            type=str,
            help="Comma-separated sample size percentages (e.g., '0.1,1,10,50')",
        )
        parser.add_argument(
            "--n-replicates",
            type=int,
            help="Number of replicates per sample size",
        )
        parser.add_argument(
            "--seed",
            type=int,
            default=42,
            help="Random seed for reproducibility (default: 42)",
        )
        parser.add_argument(
            "-o",
            "--output-file",
            help="Output CSV file (default: stdout)",
        )
        parser.add_argument(
            "--threads",
            type=int,
            default=1,
            help="Number of threads for parallel processing (default: 1)",
        )
        self._add_common_args(parser)

    def _parse_sample_sizes(self, sample_sizes_str: str) -> List[float]:
        """Parse sample size string into list of percentages.

        Parameters
        ----------
        sample_sizes_str : str
            Comma-separated sample size percentages

        Returns
        -------
        List[float]
            List of sample size percentages (0-100)
        """
        sizes = []
        for size_str in sample_sizes_str.split(","):
            size = float(size_str.strip())
            if size <= 0 or size > 100:
                raise ValueError(f"Sample size must be between 0 and 100, got {size}")
            sizes.append(size)
        return sorted(sizes)

    def _calculate_metrics_for_subsample(
        self,
        subsample_df: pd.DataFrame,
        metrics: List[Dict[str, Any]],
        sep: str,
    ) -> Dict[str, Any]:
        """Calculate all metrics for a single subsample.

        Parameters
        ----------
        subsample_df : pd.DataFrame
            Subsampled dataframe
        metrics : List[Dict[str, Any]]
            List of metric configurations
        sep : str
            CSV separator

        Returns
        -------
        Dict[str, Any]
            Dictionary of metric names to values
        """
        results = {}

        # Create temporary file for subsample
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False, newline=""
        ) as temp_file:
            temp_path = temp_file.name
            subsample_df.to_csv(temp_file, index=False, sep=sep)

        try:
            # Calculate each metric
            for metric_config in metrics:
                method = metric_config.get("method")
                name = metric_config.get("name")

                if not method or not name:
                    continue

                if method not in STAT_REGISTRY:
                    logger.warning(f"Unknown metric method '{method}', skipping")
                    continue

                stat_class = STAT_REGISTRY[method]
                step_params = {
                    k: v for k, v in metric_config.items() if k not in ["method", "name"]
                }

                try:
                    result = stat_class._from_step(temp_path, sep=sep, **step_params)
                    results[name] = result
                except Exception as e:
                    logger.error(f"Error calculating {name}: {e}")
                    results[name] = None

        finally:
            # Clean up temporary file
            if os.path.exists(temp_path):
                os.unlink(temp_path)

        return results

    def run(self) -> None:
        """Run the subsample command.

        This method orchestrates the subsampling process, handling configuration
        loading, file processing, and result output.

        Raises
        ------
        ValueError
            If required arguments are missing or invalid
        """
        logger.info("Starting subsample process")

        # Load config if provided
        if not self.args.config:
            raise ValueError("Config file is required for subsample command. Use --config to specify.")

        self._load_config(self.args.config)

        # Merge config and args with defaults
        defaults = {
            "sep": ",",
            "seed": 42,
            "threads": 1,
        }
        self._merge_config_and_args(defaults)

        # Validate input file
        if not os.path.exists(self.args.input_file):
            raise ValueError(f"Input file not found: {self.args.input_file}")

        # Get metrics - try subsample.metrics first, fall back to stats.metrics
        metrics = None
        if "subsample" in self._toml_doc and "metrics" in self._toml_doc["subsample"]:
            metrics = self._toml_doc["subsample"]["metrics"]
            logger.info(f"Using {len(metrics)} metrics from [[subsample.metrics]]")
        elif "stats" in self._toml_doc and "metrics" in self._toml_doc["stats"]:
            metrics = self._toml_doc["stats"]["metrics"]
            logger.info(f"Using {len(metrics)} metrics from [[stats.metrics]]")
        else:
            raise ValueError(
                "Config file must contain [[subsample.metrics]] or [[stats.metrics]] sections"
            )

        if not metrics:
            raise ValueError("No metrics defined in config")

        # Get sample sizes
        if not self.args.sample_sizes:
            raise ValueError("--sample-sizes is required")

        sample_sizes = self._parse_sample_sizes(self.args.sample_sizes)
        logger.info(f"Sample sizes: {sample_sizes}%")

        # Get n_replicates
        if not self.args.n_replicates:
            raise ValueError("--n-replicates is required")

        n_replicates = self.args.n_replicates
        logger.info(f"Replicates per sample size: {n_replicates}")

        # Load input data
        logger.info(f"Loading input file: {self.args.input_file}")
        df = pd.read_csv(self.args.input_file, sep=self.args.sep)
        total_rows = len(df)
        logger.info(f"Loaded {total_rows} rows")

        # Prepare results list
        all_results = []

        # Create progress bar
        total_iterations = len(sample_sizes) * n_replicates
        pbar = tqdm(
            total=total_iterations,
            desc="Subsampling",
            unit="subsample",
        )

        # Perform subsampling
        for size_idx, sample_pct in enumerate(sample_sizes):
            sample_n = max(1, int(total_rows * sample_pct / 100.0))
            
            # Ensure we don't try to sample more rows than available
            if sample_n > total_rows:
                logger.warning(
                    f"Sample size {sample_n} exceeds total rows {total_rows}, "
                    f"using {total_rows} instead"
                )
                sample_n = total_rows

            for rep_idx in range(n_replicates):
                # Calculate unique seed for this subsample
                seed = self.args.seed + (size_idx * n_replicates + rep_idx)

                # Perform subsampling
                subsample_df = df.sample(n=sample_n, random_state=seed)

                # Calculate metrics
                metric_results = self._calculate_metrics_for_subsample(
                    subsample_df, metrics, self.args.sep
                )

                # Store results for each metric
                for metric_name, metric_value in metric_results.items():
                    all_results.append({
                        "sample_size_pct": sample_pct,
                        "sample_size_n": sample_n,
                        "replicate": rep_idx,
                        "metric_name": metric_name,
                        "metric_value": metric_value,
                    })

                pbar.update(1)

        pbar.close()

        if not all_results:
            raise ValueError("No results generated")

        logger.info(f"Generated {len(all_results)} result rows")

        # Write results
        fieldnames = ["sample_size_pct", "sample_size_n", "replicate", "metric_name", "metric_value"]

        if self.args.output_file:
            logger.info(f"Writing results to {self.args.output_file}")
            with open(self.args.output_file, "w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(all_results)
        else:
            writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_results)

        logger.info("Subsample process completed successfully")


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.

