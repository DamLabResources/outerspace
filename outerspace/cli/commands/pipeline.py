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
import yaml

from snakemake.api import SnakemakeApi
from snakemake.resources import DefaultResources
from snakemake.settings.types import (
    ConfigSettings,
    DAGSettings,
    DeploymentSettings,
    ExecutionSettings,
    OutputSettings,
    RemoteExecutionSettings,
    ResourceSettings,
    SchedulingSettings,
    StorageSettings,
    WorkflowSettings,
)

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

    def _parse_snakemake_config_dict(
        self, snakemake_args: List[str]
    ) -> Dict[str, Any]:
        """Parse config values from snakemake arguments.

        Parameters
        ----------
        snakemake_args : List[str]
            Additional Snakemake arguments from user

        Returns
        -------
        Dict[str, Any]
            Dictionary of config values
        """
        config = {"toml": self.args.config_file}
        
        # Parse --config arguments from snakemake_args
        i = 0
        while i < len(snakemake_args):
            if snakemake_args[i] == "--config":
                if i + 1 < len(snakemake_args):
                    # Parse key=value
                    config_arg = snakemake_args[i + 1]
                    if "=" in config_arg:
                        key, value = config_arg.split("=", 1)
                        config[key] = value
                    i += 2
                else:
                    i += 1
            else:
                i += 1
        
        return config

    def _parse_execution_cores(self, snakemake_args: List[str]) -> Optional[int]:
        """Parse cores/threads from snakemake arguments.

        Parameters
        ----------
        snakemake_args : List[str]
            Additional Snakemake arguments from user

        Returns
        -------
        Optional[int]
            Number of cores to use, or None if not specified
        """
        # Look for --cores or -c arguments
        for i, arg in enumerate(snakemake_args):
            if arg in ["--cores", "-c"] and i + 1 < len(snakemake_args):
                try:
                    return int(snakemake_args[i + 1])
                except ValueError:
                    pass
        return None

    def _parse_jobs(self, snakemake_args: List[str]) -> Optional[int]:
        """Parse jobs/nodes from snakemake arguments for remote execution.

        Parameters
        ----------
        snakemake_args : List[str]
            Additional Snakemake arguments from user

        Returns
        -------
        Optional[int]
            Number of jobs/nodes to use, or None if not specified
        """
        # Look for --jobs or -j arguments
        for i, arg in enumerate(snakemake_args):
            if arg in ["--jobs", "-j"] and i + 1 < len(snakemake_args):
                try:
                    return int(snakemake_args[i + 1])
                except ValueError:
                    pass
        return None

    def _parse_executor(self, snakemake_args: List[str]) -> str:
        """Parse executor from snakemake arguments.

        Parameters
        ----------
        snakemake_args : List[str]
            Additional Snakemake arguments from user

        Returns
        -------
        str
            Executor name (e.g., 'local', 'slurm', 'dryrun')
        """
        # Check for dry-run first (highest priority)
        if "--dry-run" in snakemake_args or "-n" in snakemake_args:
            return "dryrun"
        
        # Look for --executor argument
        for i, arg in enumerate(snakemake_args):
            if arg == "--executor" and i + 1 < len(snakemake_args):
                return snakemake_args[i + 1]
        
        # Default to local executor
        return "local"

    def _parse_profile(self, snakemake_args: List[str]) -> Optional[Path]:
        """Parse profile path from snakemake arguments.

        Parameters
        ----------
        snakemake_args : List[str]
            Additional Snakemake arguments from user

        Returns
        -------
        Optional[Path]
            Path to profile directory, or None if not specified
        """
        # Look for --profile argument
        for i, arg in enumerate(snakemake_args):
            if arg == "--profile" and i + 1 < len(snakemake_args):
                profile_path = Path(snakemake_args[i + 1])
                if profile_path.exists():
                    return profile_path
                else:
                    logger.warning(f"Profile directory not found: {profile_path}")
                    return None
        return None

    def _load_profile_config(self, profile_path: Path) -> Dict[str, Any]:
        """Load configuration from a Snakemake profile directory.

        Parameters
        ----------
        profile_path : Path
            Path to the profile directory

        Returns
        -------
        Dict[str, Any]
            Profile configuration dictionary
        """
        # Try different config file names (v8+ first, then default)
        config_files = [
            profile_path / "config.v8+.yaml",
            profile_path / "config.v8+.yml",
            profile_path / "config.yaml",
            profile_path / "config.yml",
        ]
        
        for config_file in config_files:
            if config_file.exists():
                logger.info(f"Loading profile config from: {config_file}")
                try:
                    with open(config_file, "r") as f:
                        profile_config = yaml.safe_load(f)
                        return profile_config if profile_config else {}
                except Exception as e:
                    logger.warning(f"Failed to load profile config from {config_file}: {e}")
        
        logger.warning(f"No config file found in profile directory: {profile_path}")
        return {}

    def _parse_default_resources(self, profile_config: Dict[str, Any]) -> Optional[DefaultResources]:
        """Parse default resources from profile configuration.

        Parameters
        ----------
        profile_config : Dict[str, Any]
            Profile configuration dictionary

        Returns
        -------
        Optional[DefaultResources]
            DefaultResources object if default-resources found in profile, None otherwise
        """
        # Check for both 'default-resources' and 'default_resources' keys
        default_resources_list = profile_config.get("default-resources") or profile_config.get("default_resources")
        
        if not default_resources_list:
            return None
        
        # Convert to list if it's a single string
        if isinstance(default_resources_list, str):
            default_resources_list = [default_resources_list]
        
        if not isinstance(default_resources_list, list):
            logger.warning(f"default-resources must be a list, got {type(default_resources_list)}")
            return None
        
        try:
            default_resources = DefaultResources(args=default_resources_list)
            logger.info(f"Loaded default resources from profile: {default_resources_list}")
            return default_resources
        except Exception as e:
            logger.warning(f"Failed to parse default-resources: {e}")
            return None

    def _is_dry_run(self, snakemake_args: List[str]) -> bool:
        """Check if dry-run mode is requested.

        Parameters
        ----------
        snakemake_args : List[str]
            Additional Snakemake arguments from user

        Returns
        -------
        bool
            True if dry-run is requested
        """
        return "--dry-run" in snakemake_args or "-n" in snakemake_args

    def _execute_snakemake(
        self, snakefile_path: Path, snakemake_args: List[str]
    ) -> None:
        """Execute Snakemake workflow using the new v9 API.

        Parameters
        ----------
        snakefile_path : Path
            Path to the Snakefile
        snakemake_args : List[str]
            Additional Snakemake arguments from user

        Raises
        ------
        SystemExit
            If Snakemake execution fails
        """
        logger.info(f"Running Snakemake workflow from: {snakefile_path}")
        logger.info(f"Additional arguments: {snakemake_args}")

        try:
            # Ensure .snakemake directory exists for logs and metadata
            snakemake_dir = Path(".snakemake")
            snakemake_dir.mkdir(exist_ok=True)
            (snakemake_dir / "log").mkdir(exist_ok=True)
            
            # Parse profile if provided
            profile_path = self._parse_profile(snakemake_args)
            profile_config = {}
            if profile_path:
                profile_config = self._load_profile_config(profile_path)
                logger.info(f"Loaded profile from: {profile_path}")
            
            # Parse configuration and execution settings
            # Profile settings can be overridden by command-line args
            config_dict = self._parse_snakemake_config_dict(snakemake_args)
            
            # Get cores/jobs from args first, fall back to profile
            cores = self._parse_execution_cores(snakemake_args)
            if cores is None and "cores" in profile_config:
                cores = profile_config.get("cores")
                logger.info(f"Using cores from profile: {cores}")
            
            jobs = self._parse_jobs(snakemake_args)
            if jobs is None and "jobs" in profile_config:
                jobs = profile_config.get("jobs")
                logger.info(f"Using jobs from profile: {jobs}")
            
            # Get executor from args first, fall back to profile
            executor = self._parse_executor(snakemake_args)
            if executor == "local" and "executor" in profile_config:
                # Only override if user didn't explicitly set executor
                if "--executor" not in snakemake_args:
                    executor = profile_config.get("executor", "local")
                    logger.info(f"Using executor from profile: {executor}")
            
            # Parse default resources from profile
            default_resources = self._parse_default_resources(profile_config)
            
            # Configure resource settings based on executor type
            # For local/dryrun: use cores
            # For remote executors (slurm, etc.): use nodes (jobs) and optionally cores
            resource_kwargs = {}
            if executor in ["local", "dryrun"]:
                # Local execution: cores is required
                resource_kwargs["cores"] = cores if cores is not None else 1
            else:
                # Remote execution: nodes (jobs) is primary, cores can be per-node
                resource_kwargs["nodes"] = jobs if jobs is not None else 1
                if cores is not None:
                    resource_kwargs["cores"] = cores
            
            # Add default resources if available
            if default_resources:
                resource_kwargs["default_resources"] = default_resources
            
            logger.info(f"Executor: {executor}")
            logger.info(f"Resource settings: {resource_kwargs}")
            if default_resources:
                logger.info(f"Default resources: {default_resources.args}")
            
            # Create output settings (use defaults for now)
            output_settings = OutputSettings()
            
            # Execute with the new API
            logger.info("Starting Snakemake workflow execution")
            with SnakemakeApi(output_settings) as snakemake_api:
                # Create workflow with explicit workdir (defaults to current directory)
                workflow_api = snakemake_api.workflow(
                    resource_settings=ResourceSettings(**resource_kwargs),
                    config_settings=ConfigSettings(
                        config=config_dict,
                        configfiles=[Path(self.args.snakemake_config)],
                    ),
                    storage_settings=StorageSettings(),
                    workflow_settings=WorkflowSettings(),
                    deployment_settings=DeploymentSettings(),
                    snakefile=snakefile_path,
                    workdir=Path.cwd(),  # Explicitly set working directory
                )
                
                # Create DAG
                dag_api = workflow_api.dag(
                    dag_settings=DAGSettings(),
                )
                
                # Get executor-specific settings
                # For remote executors, the plugin's get_settings method needs to be called
                executor_settings = None
                if executor not in ["local", "dryrun", "touch"]:
                    try:
                        from snakemake_interface_executor_plugins.registry import ExecutorPluginRegistry
                        from dataclasses import fields
                        
                        executor_plugin = ExecutorPluginRegistry().get_plugin(executor)
                        
                        # Create args object with all required attributes set to None/defaults
                        # The plugin will use defaults for any None values
                        class ExecutorArgs:
                            def __init__(self):
                                # Set default attributes based on the executor's settings class
                                if executor_plugin.settings_cls:
                                    for field in fields(executor_plugin.settings_cls):
                                        # Use the prefixed name (e.g., "slurm_logdir")
                                        prefixed_name = f"{executor_plugin.cli_prefix}_{field.name}"
                                        
                                        # Check if profile has this setting
                                        profile_value = profile_config.get(prefixed_name)
                                        if profile_value is not None:
                                            setattr(self, prefixed_name, profile_value)
                                            logger.debug(f"Using {prefixed_name} from profile: {profile_value}")
                                        else:
                                            # Set to None - the plugin will use its defaults
                                            setattr(self, prefixed_name, None)
                                
                                # Set sensible defaults for common executor settings
                                # These will override None values (but not profile values)
                                if executor == "slurm":
                                    # Use the .snakemake/log directory we created
                                    if not hasattr(self, 'slurm_logdir') or self.slurm_logdir is None:
                                        self.slurm_logdir = snakemake_dir / "log"
                        
                        args_obj = ExecutorArgs()
                        executor_settings = executor_plugin.get_settings(args_obj)
                        logger.info(f"Loaded executor settings for {executor}")
                    except Exception as e:
                        logger.warning(f"Could not load executor settings for {executor}: {e}")
                        logger.warning("Continuing with default settings...")
                
                # Execute workflow
                dag_api.execute_workflow(
                    executor=executor,
                    execution_settings=ExecutionSettings(),
                    remote_execution_settings=RemoteExecutionSettings(),
                    scheduling_settings=SchedulingSettings(),
                    executor_settings=executor_settings,
                )
            
            logger.info("Pipeline completed successfully")

        except Exception as e:
            logger.error(f"Pipeline execution failed: {e}")
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

        # Parse additional Snakemake arguments
        snakemake_args = self._parse_snakemake_args(self.args.snakemake_args)

        # Locate the Snakefile
        snakefile_path = self._locate_snakefile()

        # Execute Snakemake
        self._execute_snakemake(snakefile_path, snakemake_args)


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
