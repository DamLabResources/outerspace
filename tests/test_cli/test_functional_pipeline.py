"""Functional tests for the OUTERSPACE CLI pipeline workflow"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import os
import pytest
import yaml
import logging
from outerspace.cli.main import Cli
from outerspace.cli.logging_config import setup_logging

# Set up logging
logger = setup_logging(level=logging.DEBUG)

# temp_workspace is created by conftest.py


def check_directory_contents(expected_roots, extension, path):
    """Check that the directory contains the expected files"""
    logger.debug(f"Checking directory: {path}")
    logger.debug(f"Expected files: {[root + extension for root in expected_roots]}")
    for root in expected_roots:
        file_path = os.path.join(path, root + extension)
        logger.debug(f"Checking file: {file_path}")
        assert os.path.exists(file_path), f"File {root} not found in {path}"


def test_pipeline_from_samplesheet(temp_workspace):
    """Test the pipeline command processing a samplesheet"""
    logger.info(f"Testing pipeline with samplesheet in workspace: {temp_workspace}")

    # Change to temporary workspace directory
    original_dir = os.getcwd()
    os.chdir(temp_workspace)

    try:
        pipeline_args = [
            "pipeline",
            "grnaquery.toml",
            "samplesheet_config.yaml",
            '--snakemake-args= "-c 1"',
        ]
        cli = Cli(pipeline_args)
        cli.run()

        # Verify that it made them based on the names in samplesheet.csv
        expected_files = ["M1_lib", "M2_lib", "shuffle"]
        check_directory_contents(
            expected_files, ".csv", os.path.join(temp_workspace, "findseq")
        )
        check_directory_contents(
            expected_files, ".csv", os.path.join(temp_workspace, "collapse")
        )
        check_directory_contents(
            expected_files, ".csv", os.path.join(temp_workspace, "count")
        )
        check_directory_contents(
            ["merged"], ".csv", os.path.join(temp_workspace, "merge")
        )
        check_directory_contents(
            ["stats"], ".csv", os.path.join(temp_workspace, "stats")
        )

    finally:
        # Change back to original directory
        os.chdir(original_dir)


def test_pipeline_from_directory(temp_workspace):
    """Test the pipeline command processing multiple files from a directory"""
    logger.info(f"Testing pipeline with directory in workspace: {temp_workspace}")

    # Change to temporary workspace directory
    original_dir = os.getcwd()
    os.chdir(temp_workspace)

    try:
        pipeline_args = [
            "pipeline",
            "grnaquery.toml",
            "directory_config.yaml",
            '--snakemake-args= "-c 1"',
        ]
        cli = Cli(pipeline_args)
        cli.run()

        # Verify outputs (same as basic workflow)
        expected_files = ["409-4_S1_L002", "2-G1L9-M1_S9_L001", "2-G1L9-M2_S12_L001"]
        check_directory_contents(
            expected_files, ".csv", os.path.join(temp_workspace, "findseq")
        )
        check_directory_contents(
            expected_files, ".csv", os.path.join(temp_workspace, "collapse")
        )
        check_directory_contents(
            expected_files, ".csv", os.path.join(temp_workspace, "count")
        )
        check_directory_contents(
            ["merged"], ".csv", os.path.join(temp_workspace, "merge")
        )
        check_directory_contents(
            ["stats"], ".csv", os.path.join(temp_workspace, "stats")
        )
    finally:
        # Change back to original directory
        os.chdir(original_dir)


def test_pipeline_from_directory_accepts_snakemake_args(temp_workspace):
    """Test the pipeline command processing multiple files from a directory"""
    logger.info(f"Testing pipeline with directory in workspace: {temp_workspace}")

    # Change to temporary workspace directory
    original_dir = os.getcwd()
    os.chdir(temp_workspace)

    try:
        pipeline_args = [
            "pipeline",
            "grnaquery.toml",
            "directory_config.yaml",
            '--snakemake-args="--dry-run"',
        ]
        cli = Cli(pipeline_args)
        cli.run()

        # Verify that the pipeline didn't run by making sure no files are in the directories
        assert not os.listdir("findseq"), "findseq should not have any files"
        assert not os.listdir("collapse"), "collapse should not have any files"
        assert not os.listdir("count"), "count should not have any files"

    finally:
        # Change back to original directory
        os.chdir(original_dir)

# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.