"""Tests for the stats command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import os
import pytest
import tempfile
import csv
from outerspace.cli.main import Cli


def test_stats_initialization():
    """Test that stats command initializes correctly with config"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file with new stepwise format
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[[stats.metrics]]
method = "simpson_diversity"
key_column = "protospacer"
barcode_column = "UMI_5prime_UMI_3prime_corrected_count"
name = "simpson"

[[stats.metrics]]
method = "shannon_diversity"
key_column = "protospacer"
barcode_column = "UMI_5prime_UMI_3prime_corrected_count"
name = "shannon"
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,UMI_5prime_UMI_3prime_corrected_count\n")
            f.write("A,1\n")
            f.write("B,2\n")
            f.write("C,3\n")

        args = ["stats", "--config", config_file, input_file]
        cli = Cli(args)
        cli.run()


def test_stats_with_multiple_metrics():
    """Test stats command with multiple different metrics"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file with multiple metrics
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[[stats.metrics]]
method = "gini_coefficient"
key_column = "key"
barcode_column = "count"
name = "gini"

[[stats.metrics]]
method = "shannon_diversity"
key_column = "key"
barcode_column = "count"
name = "shannon"
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("key,count\n")
            f.write("A,10\n")
            f.write("B,5\n")
            f.write("C,3\n")

        args = ["stats", "--config", config_file, input_file]
        cli = Cli(args)
        cli.run()


def test_stats_missing_required_args():
    """Test that stats command requires config file"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("key,count\n")
            f.write("A,1\n")

        args = ["stats", input_file]
        with pytest.raises(
            ValueError, match="Config file is required for stats command"
        ):
            cli = Cli(args)
            cli.run()


def test_stats_missing_metrics():
    """Test that stats command requires metrics in config"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file without stats.metrics
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[count]
key_column = "protospacer"
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,count\n")
            f.write("A,1\n")

        args = ["stats", "--config", config_file, input_file]
        with pytest.raises(ValueError, match="Config file must contain \\[\\[stats\\.metrics\\]\\] sections"):
            cli = Cli(args)
            cli.run()


def test_stats_nonexistent_input_file():
    """Test that stats command handles nonexistent input file"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[[stats.metrics]]
method = "gini_coefficient"
key_column = "protospacer"
barcode_column = "count"
name = "gini"
"""
            )

        args = ["stats", "--config", config_file, "nonexistent.csv"]
        with pytest.raises(ValueError, match="No input files found"):
            cli = Cli(args)
            cli.run()


def test_stats_nonexistent_columns():
    """Test that stats command handles nonexistent columns"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file with nonexistent columns
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[[stats.metrics]]
method = "gini_coefficient"
key_column = "nonexistent"
barcode_column = "also_nonexistent"
name = "gini"
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,UMI_5prime_UMI_3prime_corrected_count\n")
            f.write("A,1\n")

        args = [
            "stats",
            "--config",
            config_file,
            input_file,
        ]
        with pytest.raises(ValueError):
            cli = Cli(args)
            cli.run()


def test_stats_with_allowed_list():
    """Test that stats command works with allowed list in metric config"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create allowed list file
        allowed_list = os.path.join(temp_dir, "allowed.txt")
        with open(allowed_list, "w") as f:
            f.write("A\nB\n")

        # Create config file with metric using allowed_list
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                f"""[[stats.metrics]]
method = "umi_recovery_rate"
key_column = "protospacer"
barcode_column = "count"
allowed_list = "{allowed_list}"
name = "recovery_rate"
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,count\n")
            f.write("A,1\n")
            f.write("B,2\n")
            f.write("C,3\n")

        args = [
            "stats",
            "--config",
            config_file,
            input_file,
        ]
        cli = Cli(args)
        cli.run()


def test_stats_multiple_input_files():
    """Test that stats command handles multiple input files"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[[stats.metrics]]
method = "simpson_diversity"
key_column = "protospacer"
barcode_column = "count"
name = "simpson"
"""
            )

        # Create test input files
        input_files = []
        for i in range(3):
            input_file = os.path.join(temp_dir, f"test{i}.csv")
            with open(input_file, "w") as f:
                f.write("protospacer,count\n")
                f.write(f"A{i},1\n")
                f.write(f"B{i},2\n")
            input_files.append(input_file)

        args = ["stats", "--config", config_file] + input_files
        cli = Cli(args)
        cli.run()


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
