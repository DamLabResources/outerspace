"""Statistics calculation command for UMI analysis.

This module provides the StatsCommand class for calculating comprehensive statistics
from UMI count data using a config-driven, stepwise approach. Each metric is
defined separately in the configuration with its own parameters.
"""

import csv
import logging
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Type
from argparse import ArgumentParser

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
STAT_REGISTRY: Dict[str, Type[BaseStatistic]] = {
    "gini_coefficient": GiniCoefficient,
    "shannon_diversity": ShannonDiversity,
    "simpson_diversity": SimpsonDiversity,
    "hill_number": HillNumber,
    "umi_recovery_rate": UMIRecoveryRate,
    "umi_efficiency_rate": UMIEfficiencyRate,
    "error_rate": ErrorRate,
}


class StatsCommand(BaseCommand):
    """Command for calculating statistics using config-defined metrics.

    This command computes statistics for UMI data using a stepwise, config-driven
    approach where each metric is defined separately with its own parameters.
    Metrics are specified in the config file under [[stats.metrics]] sections.
    """

    def _init_parser(self, subparser: ArgumentParser) -> None:
        """Initialize command-specific argument parser.

        Parameters
        ----------
        subparser : ArgumentParser
            Subparser group to add command arguments to
        """
        parser = subparser.add_parser(
            "stats",
            help="Calculate statistics using config-defined metrics",
        )
        parser.add_argument(
            "input_files",
            nargs="+",
            help="Input CSV file(s) to process (supports glob patterns)",
        )
        parser.add_argument("--sep", default=",", help="CSV separator (default: ,)")
        parser.add_argument(
            "--threads",
            type=int,
            default=1,
            help="Number of threads for parallel file processing (default: 1)",
        )
        parser.add_argument(
            "-o",
            "--output-file",
            help="Output CSV file (default: stdout)",
        )
        self._add_common_args(parser)

    def _calculate_stats_for_file(
        self,
        input_file: str,
        metrics: List[Dict[str, Any]],
        sep: str = ",",
        pbar: Optional[tqdm] = None,
    ) -> Dict[str, Any]:
        """Calculate all configured metrics for a single input file.

        Parameters
        ----------
        input_file : str
            Path to input CSV file
        metrics : List[Dict[str, Any]]
            List of metric configurations from config file
        sep : str, default=','
            CSV separator
        pbar : Optional[tqdm], default=None
            Optional progress bar to update after processing

        Returns
        -------
        Dict[str, Any]
            Dictionary containing all calculated statistics

        Raises
        ------
        ValueError
            If metric method is not found in registry or calculation fails
        """
        filename = os.path.basename(input_file)
        
        # Use tqdm.write for logging to avoid conflicts with progress bars
        if pbar is not None:
            tqdm.write(f"Processing {filename}...")
        else:
            logger.info(f"Calculating statistics for {filename}")

        stats = {"filename": filename}

        # Create a nested progress bar for metrics if we're showing progress
        metric_pbar = None
        if pbar is not None and len(metrics) > 1:
            metric_pbar = tqdm(
                total=len(metrics),
                desc=f"  Metrics",
                leave=False,
                position=1,
            )

        for metric_config in metrics:
            method = metric_config.get("method")
            name = metric_config.get("name")

            if not method:
                raise ValueError("Metric configuration must include 'method' field")
            if not name:
                raise ValueError("Metric configuration must include 'name' field")

            if method not in STAT_REGISTRY:
                raise ValueError(
                    f"Unknown metric method '{method}'. "
                    f"Available methods: {', '.join(STAT_REGISTRY.keys())}"
                )

            # Get the statistic class
            stat_class = STAT_REGISTRY[method]

            # Extract parameters for this metric (excluding method and name)
            step_params = {k: v for k, v in metric_config.items() if k not in ["method", "name"]}

            try:
                # Calculate the statistic using _from_step
                logger.debug(f"Calling {method} with params: {step_params}")
                result = stat_class._from_step(input_file, sep=sep, **step_params)
                stats[name] = result

                # Log result
                if pbar is not None:
                    # Use tqdm.write to avoid conflicts with progress bars
                    if result is not None:
                        if isinstance(result, dict):
                            tqdm.write(f"    {name}: {result}")
                        else:
                            tqdm.write(f"    {name}: {result:.6f}")
                    else:
                        tqdm.write(f"    {name}: None")
                else:
                    if result is not None:
                        if isinstance(result, dict):
                            logger.info(f"  {name}: {result}")
                        else:
                            logger.info(f"  {name}: {result:.6f}")
                    else:
                        logger.info(f"  {name}: None")
                
                # Update metric progress bar
                if metric_pbar is not None:
                    metric_pbar.update(1)

            except Exception as e:
                if pbar is not None:
                    tqdm.write(f"ERROR: Failed to calculate {name}: {e}")
                logger.error(f"Error calculating {name} ({method}) with params {step_params}: {e}")
                if metric_pbar is not None:
                    metric_pbar.close()
                raise ValueError(f"Failed to calculate {name}: {e}")

        # Close metric progress bar if it was created
        if metric_pbar is not None:
            metric_pbar.close()
        
        # Update main progress bar
        if pbar is not None:
            pbar.update(1)

        return stats

    def run(self) -> None:
        """Run the stats command.

        This method orchestrates the statistics calculation process, handling
        configuration loading, file processing, and result output with comprehensive
        error handling and logging.

        Raises
        ------
        ValueError
            If required arguments are missing or invalid, or if no files are
            successfully processed
        """
        logger.info("Starting statistics calculation process")

        # Load config if provided
        if not self.args.config:
            raise ValueError("Config file is required for stats command. Use --config to specify.")

        self._load_config(self.args.config)

        # Merge config and args with defaults
        defaults = {"sep": ","}
        self._merge_config_and_args(defaults)

        # Get metrics from config
        # Access the full TOML document to get stats.metrics
        if "stats" not in self._toml_doc or "metrics" not in self._toml_doc["stats"]:
            raise ValueError(
                "Config file must contain [[stats.metrics]] sections defining the metrics to calculate"
            )

        metrics = self._toml_doc["stats"]["metrics"]
        if not metrics:
            raise ValueError("No metrics defined in [[stats.metrics]] sections")

        logger.info(f"Found {len(metrics)} metrics to calculate")

        try:
            # Filter to existing files
            input_files = [f for f in self.args.input_files if os.path.exists(f)]
            
            # Warn about missing files
            for f in self.args.input_files:
                if not os.path.exists(f):
                    logger.warning(f"Input file not found: {f}")
            
            if not input_files:
                raise ValueError("No input files found")
            
            # Process files with optional multi-threading
            all_stats = []
            
            # Create progress bar for file processing
            file_pbar = tqdm(
                total=len(input_files),
                desc="Processing files",
                unit="file",
                position=0,
            )
            
            if self.args.threads > 1:
                logger.info(f"Processing {len(input_files)} files with {self.args.threads} threads")
                
                # Use ThreadPoolExecutor for parallel processing
                with ThreadPoolExecutor(max_workers=self.args.threads) as executor:
                    # Submit all tasks
                    future_to_file = {
                        executor.submit(
                            self._calculate_stats_for_file,
                            input_file,
                            metrics,
                            self.args.sep,
                            file_pbar
                        ): input_file
                        for input_file in input_files
                    }
                    
                    # Collect results as they complete
                    for future in as_completed(future_to_file):
                        input_file = future_to_file[future]
                        try:
                            stats = future.result()
                            all_stats.append(stats)
                        except Exception as e:
                            tqdm.write(f"ERROR: Error processing {input_file}: {e}")
                            logger.error(f"Error processing {input_file}: {e}")
            else:
                # Single-threaded processing
                for input_file in input_files:
                    try:
                        stats = self._calculate_stats_for_file(
                            input_file,
                            metrics,
                            sep=self.args.sep,
                            pbar=file_pbar,
                        )
                        all_stats.append(stats)
                    except Exception as e:
                        tqdm.write(f"ERROR: Error processing {input_file}: {e}")
                        logger.error(f"Error processing {input_file}: {e}")
                        file_pbar.update(1)
                        continue
            
            # Close the progress bar
            file_pbar.close()

            if not all_stats:
                raise ValueError("No files were successfully processed")

            logger.info(f"Successfully processed {len(all_stats)} files")

            # Write all results to output file or stdout
            if self.args.output_file:
                logger.info(f"Writing results to {self.args.output_file}")
                with open(self.args.output_file, 'w', newline='') as f:
                    writer = csv.DictWriter(f, fieldnames=all_stats[0].keys())
                    writer.writeheader()
                    writer.writerows(all_stats)
            else:
                writer = csv.DictWriter(sys.stdout, fieldnames=all_stats[0].keys())
                writer.writeheader()
                writer.writerows(all_stats)

            logger.info("Statistics calculation completed successfully")

        except Exception as e:
            logger.error(f"Error calculating statistics: {e}")
            raise ValueError(f"Error calculating statistics: {e}")


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
