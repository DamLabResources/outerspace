"""Pipeline command for running complete OUTERSPACE workflows.

This module provides the PipelineCommand class for executing complete OUTERSPACE
pipelines using Snakemake. It handles configuration loading, argument parsing,
and workflow execution with comprehensive error handling.
"""

import logging
import os
import shlex
import sys
from pathlib import Path
try:
    # Python 3.9+
    from importlib.resources import files as resource_files, as_file  # type: ignore
    _HAS_FILES_API = True
except Exception:  # pragma: no cover - fallback for Python 3.8
    import importlib.resources as _resources  # type: ignore
    resource_files = None  # type: ignore
    as_file = None  # type: ignore
    _HAS_FILES_API = False
from typing import Any, Dict, List, Optional
from argparse import ArgumentParser

import pulp
import snakemake
import yaml

from outerspace.cli.commands.base import BaseCommand

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

# Monkey patch pulp.list_solvers if it doesn't exist
if not hasattr(pulp, "list_solvers"):
    pulp.list_solvers = pulp.listSolvers


class PipelineCommand(BaseCommand):
    """Command for running the complete OUTERSPACE pipeline using Snakemake.

    This command orchestrates the execution of complete OUTERSPACE workflows
    using Snakemake as the workflow engine. It handles configuration management,
    argument parsing, and workflow execution with comprehensive error handling.
    """

    def __init__(self, args: Optional[Any] = None) -> None:
        """Initialize the pipeline command.

        Parameters
        ----------
        args : Optional[Any], default=None
            Parsed command-line arguments
        """
        super().__init__(args=args)

    def _init_parser(self, subparser: ArgumentParser) -> None:
        """Initialize command-specific argument parser.

        Parameters
        ----------
        subparser : ArgumentParser
            Subparser group to add command arguments to
        """
        parser = subparser.add_parser(
            "pipeline", help="Run the complete OUTERSPACE pipeline using Snakemake"
        )

        # Required arguments
        parser.add_argument(
            "config_file", help="TOML configuration file with search patterns"
        )
        parser.add_argument(
            "snakemake_config", help="YAML configuration file for Snakemake workflow"
        )

        # Optional arguments
        parser.add_argument(
            "--snakemake-args",
            help='Additional arguments to pass to Snakemake (e.g. --snakemake-args="--dry-run --cores 4")',
        )

    def _load_snakemake_config(self, config_file: str) -> Dict[str, Any]:
        """Load Snakemake configuration from YAML file.

        Parameters
        ----------
        config_file : str
            Path to Snakemake configuration file

        Returns
        -------
        Dict[str, Any]
            Loaded configuration dictionary

        Raises
        ------
        ValueError
            If configuration file doesn't exist or is invalid
        """
        logger.info(f"Loading Snakemake config from: {config_file}")

        if not os.path.exists(config_file):
            raise ValueError(f"Snakemake configuration file not found: {config_file}")

        try:
            with open(config_file, "r") as f:
                config = yaml.safe_load(f)
                logger.debug(f"Loaded config: {config}")
                return config
        except Exception as e:
            logger.error(f"Failed to load Snakemake config from {config_file}: {e}")
            raise ValueError(f"Invalid Snakemake configuration file {config_file}: {e}")

    def _parse_snakemake_args(self, args_str: Optional[str]) -> List[str]:
        """Parse additional Snakemake arguments.

        Parameters
        ----------
        args_str : Optional[str]
            String of arguments to parse (e.g. "--dry-run --cores 4")

        Returns
        -------
        List[str]
            List of parsed arguments

        Notes
        -----
        Uses shlex.split to handle quoted arguments correctly and preserve
        argument boundaries.
        """
        if not args_str:
            return []

        # Parse the arguments using shlex to handle quotes correctly
        parsed_args = shlex.split(args_str)
        logger.debug(f"Parsed Snakemake arguments: {parsed_args}")
        return parsed_args

    def _get_packaged_snakefile(self) -> Optional[Path]:
        """Attempt to locate the packaged Snakefile from installed package.

        Returns
        -------
        Optional[Path]
            Path to the packaged Snakefile if found, None otherwise

        Notes
        -----
        Tries Python 3.9+ importlib.resources.files API first, then falls back
        to Python 3.8 compatible method using package __file__ location.
        """
        try:
            if _HAS_FILES_API:
                # Python 3.9+ method
                resource = resource_files("outerspace").joinpath("workflow", "Snakefile")
                # Test if resource exists
                if resource.is_file():
                    return Path(str(resource))
            else:
                # Python 3.8 fallback: resolve path under installed outerspace package directory
                import outerspace  # type: ignore
                snakefile_path = Path(outerspace.__file__).resolve().parent / "workflow" / "Snakefile"
                if snakefile_path.exists():
                    return snakefile_path
        except Exception as e:
            logger.debug(f"Failed to locate packaged Snakefile: {e}")
        
        return None

    def _get_repo_snakefile(self) -> Optional[Path]:
        """Attempt to locate the Snakefile in the repository directory.

        Returns
        -------
        Optional[Path]
            Path to the repository Snakefile if found, None otherwise

        Notes
        -----
        This is used as a fallback during development when the package
        is not installed or when running from source.
        """
        # Go up from cli/commands/pipeline.py -> cli/commands -> cli -> outerspace -> repo root
        repo_snakefile = Path(__file__).resolve().parents[3] / "workflow" / "Snakefile"
        if repo_snakefile.exists():
            logger.debug(f"Found repository Snakefile at: {repo_snakefile}")
            return repo_snakefile
        return None

    def _locate_snakefile(self) -> Path:
        """Locate the Snakefile to use for the pipeline.

        Returns
        -------
        Path
            Path to the Snakefile

        Raises
        ------
        ValueError
            If Snakefile cannot be found in package or repository

        Notes
        -----
        Attempts to locate Snakefile in the following order:
        1. Installed package location (for pip-installed package)
        2. Repository location (for development/editable install)
        """
        # Try packaged Snakefile first
        snakefile_path = self._get_packaged_snakefile()
        if snakefile_path:
            logger.info(f"Using packaged Snakefile: {snakefile_path}")
            return snakefile_path

        # Fall back to repository Snakefile
        snakefile_path = self._get_repo_snakefile()
        if snakefile_path:
            logger.info(f"Using repository Snakefile: {snakefile_path}")
            return snakefile_path

        # Unable to locate Snakefile
        raise ValueError(
            "Unable to locate Snakefile. Please ensure the package is properly "
            "installed or you are running from the repository root."
        )

    def _build_snakemake_argv(
        self, snakefile_path: Path, snakemake_args: List[str]
    ) -> List[str]:
        """Build the argument vector for Snakemake execution.

        Parameters
        ----------
        snakefile_path : Path
            Path to the Snakefile
        snakemake_args : List[str]
            Additional Snakemake arguments from user

        Returns
        -------
        List[str]
            Complete argument list for Snakemake
        """
        # Start with the program name and Snakefile
        new_argv = ["snakemake", "-s", str(snakefile_path)]

        # Add config file if not already specified
        if "--configfile" not in snakemake_args:
            new_argv.extend(["--configfile", self.args.snakemake_config])

        # Add remaining arguments
        new_argv.extend(snakemake_args)

        # Add our config as a config value
        if "--config" not in snakemake_args:
            new_argv.extend(["--config", f"toml={self.args.config_file}"])

        return new_argv

    def _execute_snakemake(self, argv: List[str]) -> None:
        """Execute Snakemake with the given arguments.

        Parameters
        ----------
        argv : List[str]
            Argument vector for Snakemake (including program name)

        Raises
        ------
        SystemExit
            If Snakemake execution fails (exit code != 0)
        """
        logger.info(f"Running Snakemake with args: {argv}")

        try:
            logger.info("Starting Snakemake workflow execution")
            snakemake.main(argv[1:])  # Skip the 'snakemake' program name
            logger.info("Pipeline completed successfully")

        except SystemExit as e:
            if e.code != 0:
                logger.error(f"Pipeline failed with exit code {e.code}")
                sys.exit(1)
            else:
                logger.info("Pipeline completed successfully")
        except Exception as e:
            logger.error(f"Unexpected error during pipeline execution: {e}")
            raise

    def run(self) -> None:
        """Run the pipeline command using Snakemake.

        This method orchestrates the pipeline execution by loading configurations,
        preparing Snakemake arguments, and executing the workflow with comprehensive
        error handling and logging.

        Raises
        ------
        ValueError
            If required configuration files are missing or invalid
        SystemExit
            If Snakemake execution fails (exit code != 0)
        """
        logger.info(f"Running pipeline with config file: {self.args.config_file}")
        logger.info(f"Snakemake config file: {self.args.snakemake_config}")

        # Load TOML config
        self._load_config(self.args.config_file)

        # Load Snakemake config
        snakemake_config = self._load_snakemake_config(self.args.snakemake_config)

        # Prepare Snakemake configuration
        config = {"toml": self.args.config_file, **snakemake_config}

        # Parse additional Snakemake arguments
        snakemake_args = self._parse_snakemake_args(self.args.snakemake_args)

        # Locate the Snakefile
        snakefile_path = self._locate_snakefile()

        # Build Snakemake arguments
        argv = self._build_snakemake_argv(snakefile_path, snakemake_args)

        # Execute Snakemake
        self._execute_snakemake(argv)


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
