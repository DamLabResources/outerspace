"""Tests for the subsample command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import os
import pytest
import tempfile
import csv
from outerspace.cli.main import Cli


def test_subsample_initialization():
    """Test that subsample command initializes correctly with config"""
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

        # Create test input file with enough rows
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,count\n")
            for i in range(100):
                f.write(f"SEQ{i},1\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "subsample",
            "--config",
            config_file,
            "--sample-sizes",
            "10,50",
            "--n-replicates",
            "2",
            "-o",
            output_file,
            input_file,
        ]
        cli = Cli(args)
        cli.run()

        # Verify output file exists
        assert os.path.exists(output_file)


def test_subsample_with_multiple_sample_sizes():
    """Test subsample command with multiple sample sizes"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[[stats.metrics]]
method = "gini_coefficient"
key_column = "key"
barcode_column = "count"
name = "gini"
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("key,count\n")
            for i in range(100):
                f.write(f"K{i},{i+1}\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "subsample",
            "--config",
            config_file,
            "--sample-sizes",
            "1,10,50",
            "--n-replicates",
            "2",
            "-o",
            output_file,
            input_file,
        ]
        cli = Cli(args)
        cli.run()

        # Verify output has correct number of rows (3 sizes x 2 replicates x 1 metric = 6)
        with open(output_file, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            assert len(rows) == 6


def test_subsample_with_multiple_metrics():
    """Test subsample command with multiple metrics"""
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
            for i in range(100):
                f.write(f"K{i},{i+1}\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "subsample",
            "--config",
            config_file,
            "--sample-sizes",
            "10,50",
            "--n-replicates",
            "2",
            "-o",
            output_file,
            input_file,
        ]
        cli = Cli(args)
        cli.run()

        # Verify output has correct number of rows (2 sizes x 2 replicates x 2 metrics = 8)
        with open(output_file, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            assert len(rows) == 8


def test_subsample_seeding_reproducibility():
    """Test that subsample command produces reproducible results with same seed"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[[stats.metrics]]
method = "gini_coefficient"
key_column = "key"
barcode_column = "count"
name = "gini"
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("key,count\n")
            for i in range(100):
                f.write(f"K{i},{i+1}\n")

        # Run twice with same seed
        output_file1 = os.path.join(temp_dir, "output1.csv")
        output_file2 = os.path.join(temp_dir, "output2.csv")

        for output_file in [output_file1, output_file2]:
            args = [
                "subsample",
                "--config",
                config_file,
                "--sample-sizes",
                "10,50",
                "--n-replicates",
                "2",
                "--seed",
                "12345",
                "-o",
                output_file,
                input_file,
            ]
            cli = Cli(args)
            cli.run()

        # Compare outputs
        with open(output_file1, "r") as f1, open(output_file2, "r") as f2:
            rows1 = list(csv.DictReader(f1))
            rows2 = list(csv.DictReader(f2))

            assert len(rows1) == len(rows2)
            for r1, r2 in zip(rows1, rows2):
                assert r1 == r2


def test_subsample_missing_config():
    """Test that subsample command requires config file"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("key,count\n")
            f.write("A,1\n")

        args = [
            "subsample",
            "--sample-sizes",
            "10",
            "--n-replicates",
            "2",
            input_file,
        ]
        with pytest.raises(
            ValueError, match="Config file is required for subsample command"
        ):
            cli = Cli(args)
            cli.run()


def test_subsample_missing_metrics():
    """Test that subsample command falls back to stats.metrics"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file with stats.metrics (no subsample.metrics)
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[[stats.metrics]]
method = "gini_coefficient"
key_column = "key"
barcode_column = "count"
name = "gini"
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("key,count\n")
            for i in range(50):
                f.write(f"K{i},{i+1}\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "subsample",
            "--config",
            config_file,
            "--sample-sizes",
            "10",
            "--n-replicates",
            "2",
            "-o",
            output_file,
            input_file,
        ]
        cli = Cli(args)
        cli.run()

        # Should succeed using stats.metrics
        assert os.path.exists(output_file)


def test_subsample_output_format():
    """Test that subsample command produces correct output format"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[[stats.metrics]]
method = "simpson_diversity"
key_column = "key"
barcode_column = "count"
name = "simpson"
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("key,count\n")
            for i in range(50):
                f.write(f"K{i},{i+1}\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "subsample",
            "--config",
            config_file,
            "--sample-sizes",
            "10,50",
            "--n-replicates",
            "2",
            "-o",
            output_file,
            input_file,
        ]
        cli = Cli(args)
        cli.run()

        # Verify output format
        with open(output_file, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)

            # Check columns
            assert reader.fieldnames == [
                "sample_size_pct",
                "sample_size_n",
                "replicate",
                "metric_name",
                "metric_value",
            ]

            # Check data
            assert len(rows) == 4  # 2 sizes x 2 replicates x 1 metric

            # Check sample sizes
            sample_sizes_pct = {float(row["sample_size_pct"]) for row in rows}
            assert sample_sizes_pct == {10.0, 50.0}

            # Check replicates
            for size_pct in [10.0, 50.0]:
                size_rows = [r for r in rows if float(r["sample_size_pct"]) == size_pct]
                replicates = {int(r["replicate"]) for r in size_rows}
                assert replicates == {0, 1}

            # Check metric name
            assert all(row["metric_name"] == "simpson" for row in rows)


def test_subsample_nonexistent_input():
    """Test that subsample command handles nonexistent input file"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[[stats.metrics]]
method = "gini_coefficient"
key_column = "key"
barcode_column = "count"
name = "gini"
"""
            )

        args = [
            "subsample",
            "--config",
            config_file,
            "--sample-sizes",
            "10",
            "--n-replicates",
            "2",
            "nonexistent.csv",
        ]
        with pytest.raises(ValueError, match="Input file not found"):
            cli = Cli(args)
            cli.run()


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.

