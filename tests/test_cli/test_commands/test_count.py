"""Tests for the count command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
import os
import tempfile
from outerspace.cli.main import Cli


def test_count_initialization():
    """Test that count command initializes correctly"""
    args = [
        "count",
        "--input-dir",
        "test_input",
        "--output-dir",
        "test_output",
        "--barcode-column",
        "umi3_umi5_corrected",
        "--key-column",
        "protospacer",
    ]
    cli = Cli(args)
    assert cli.args.input_dir == "test_input"
    assert cli.args.output_dir == "test_output"
    assert cli.args.barcode_column == "umi3_umi5_corrected"
    assert cli.args.key_column == "protospacer"
    assert cli.args.sep == ","
    assert cli.args.detailed is False


def test_count_single_file_initialization():
    """Test that count command initializes correctly in single file mode"""
    args = [
        "count",
        "--input-file",
        "test_input.csv",
        "--output-file",
        "test_output.csv",
        "--barcode-column",
        "umi3_umi5_corrected",
        "--key-column",
        "protospacer",
    ]
    cli = Cli(args)
    assert cli.args.input_file == "test_input.csv"
    assert cli.args.output_file == "test_output.csv"
    assert cli.args.barcode_column == "umi3_umi5_corrected"
    assert cli.args.key_column == "protospacer"
    assert cli.args.sep == ","
    assert cli.args.detailed is False


def test_count_missing_input():
    """Test that count command handles missing input"""
    args = [
        "count",
        "--output-dir",
        "test_output",
        "--barcode-column",
        "umi3_umi5_corrected",
        "--key-column",
        "protospacer",
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2


def test_count_missing_output():
    """Test that count command handles missing output"""
    args = [
        "count",
        "--input-dir",
        "test_input",
        "--barcode-column",
        "umi3_umi5_corrected",
        "--key-column",
        "protospacer",
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2


def test_count_missing_barcode_column():
    """Test that count command handles missing barcode column"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[count]
barcode_column = "UMI_5prime_UMI_3prime_corrected"
key_column = "protospacer"
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("header1,header2\nvalue1,value2\n")

        args = [
            "count",
            "--input-file",
            input_file,
            "--output-file",
            os.path.join(temp_dir, "output.csv"),
            "--key-column",
            "protospacer",
        ]
        with pytest.raises(
            ValueError, match="Please provide either --barcode-column or --config"
        ):
            cli = Cli(args)
            cli.run()


def test_count_missing_key_column():
    """Test that count command handles missing key column"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[count]
barcode_column = "UMI_5prime_UMI_3prime_corrected"
key_column = "protospacer"
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("header1,header2\nvalue1,value2\n")

        args = [
            "count",
            "--input-file",
            input_file,
            "--output-file",
            os.path.join(temp_dir, "output.csv"),
            "--barcode-column",
            "umi3_umi5_corrected",
        ]
        with pytest.raises(
            ValueError, match="Please provide either --key-column or --config"
        ):
            cli = Cli(args)
            cli.run()


def test_count_invalid_downsample():
    """Test that count command handles invalid downsample value"""
    args = [
        "count",
        "--input-dir",
        "test_input",
        "--output-dir",
        "test_output",
        "--barcode-column",
        "umi3_umi5_corrected",
        "--key-column",
        "protospacer",
        "--downsample",
        "-1",
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()


def test_count_nonexistent_input_dir():
    """Test that count command handles nonexistent input directory"""
    args = [
        "count",
        "--input-dir",
        "nonexistent_dir",
        "--output-dir",
        "test_output",
        "--barcode-column",
        "umi3_umi5_corrected",
        "--key-column",
        "protospacer",
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()


def test_count_empty_input_dir():
    """Test that count command handles empty input directory"""
    with tempfile.TemporaryDirectory() as temp_dir:
        args = [
            "count",
            "--input-dir",
            temp_dir,
            "--output-dir",
            "test_output",
            "--barcode-column",
            "umi3_umi5_corrected",
            "--key-column",
            "protospacer",
        ]
        cli = Cli(args)
        with pytest.raises(ValueError):
            cli.run()


def test_count_single_file_missing_output_file():
    """Test that count command handles missing output file for single file mode"""
    args = [
        "count",
        "--input-file",
        "test_input.csv",
        "--barcode-column",
        "umi3_umi5_corrected",
        "--key-column",
        "protospacer",
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2


def test_count_single_file_nonexistent_input():
    """Test that count command handles nonexistent input file in single file mode"""
    args = [
        "count",
        "--input-file",
        "nonexistent.csv",
        "--output-file",
        "test_output.csv",
        "--barcode-column",
        "umi3_umi5_corrected",
        "--key-column",
        "protospacer",
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()


def test_count_single_file_basic():
    """Test basic single file counting functionality"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,umi3_umi5_corrected\n")
            f.write("key1,bar1\n")
            f.write("key1,bar2\n")
            f.write("key2,bar3\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "count",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--barcode-column",
            "umi3_umi5_corrected",
            "--key-column",
            "protospacer",
        ]
        cli = Cli(args)
        cli.run()

        # Verify output file exists and has correct content
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            content = f.read()
            assert "key1,2" in content
            assert "key2,1" in content


def test_count_single_file_with_allowed_list():
    """Test single file counting with allowed list"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,umi3_umi5_corrected\n")
            f.write("key1,bar1\n")
            f.write("key1,bar2\n")
            f.write("key2,bar3\n")
            f.write("key3,bar4\n")  # This key should be filtered out

        # Create allowed list file
        allowed_list = os.path.join(temp_dir, "allowed.txt")
        with open(allowed_list, "w") as f:
            f.write("key1\n")
            f.write("key2\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "count",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--barcode-column",
            "umi3_umi5_corrected",
            "--key-column",
            "protospacer",
            "--allowed-list",
            allowed_list,
        ]
        cli = Cli(args)
        cli.run()

        # Verify output file exists and has correct content
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            content = f.read()
            assert "key1,2" in content
            assert "key2,1" in content
            assert "key3" not in content


def test_count_single_file_with_allowed_list_no_rescue_when_disabled():
    """Near-miss keys should not be counted if rescue is disabled"""
    with tempfile.TemporaryDirectory() as temp_dir:
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,umi3_umi5_corrected\n")
            f.write("key1,bar1\n")
            f.write("keyX,bar2\n")  # Near 'key1'

        allowed_list = os.path.join(temp_dir, "allowed.txt")
        with open(allowed_list, "w") as f:
            f.write("key1\n")

        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "count",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--barcode-column",
            "umi3_umi5_corrected",
            "--key-column",
            "protospacer",
            "--allowed-list",
            allowed_list,
            # rescue NOT enabled
        ]
        cli = Cli(args)
        cli.run()

        with open(output_file, "r") as f:
            content = f.read()
            assert "key1,1" in content
            assert "keyX" not in content




def test_count_with_config_file():
    """Test that count command loads and uses config file correctly"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[count]
sep = ";"
detailed = true
downsample = 0.5
random_seed = 42
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer;umi3_umi5_corrected\n")
            f.write("key1;bar1\n")
            f.write("key1;bar2\n")
            f.write("key2;bar3\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "count",
            "--config",
            config_file,
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--barcode-column",
            "umi3_umi5_corrected",
            "--key-column",
            "protospacer",
        ]
        cli = Cli(args)
        cli.run()

        # Verify config values were applied
        assert cli.args.sep == ";"
        assert cli.args.detailed is True
        assert cli.args.downsample == 0.5
        assert cli.args.random_seed == 42


def test_count_config_override():
    """Test that command line arguments override config file settings"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[count]
sep = ";"
detailed = true
downsample = 0.5
random_seed = 42
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer|umi3_umi5_corrected\n")
            f.write("key1|bar1\n")
            f.write("key1|bar2\n")
            f.write("key2|bar3\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "count",
            "--config",
            config_file,
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--barcode-column",
            "umi3_umi5_corrected",
            "--key-column",
            "protospacer",
            "--sep",
            "|",  # Override config
            "--downsample",
            "0.25",  # Override config
        ]
        cli = Cli(args)
        cli.run()

        # Verify command line args took precedence
        assert cli.args.sep == "|"
        assert cli.args.downsample == 0.25
        assert cli.args.detailed is True  # This one should come from config
        assert cli.args.random_seed == 42  # This one should come from config


def test_count_nonexistent_config():
    """Test that count command handles nonexistent config file"""
    args = [
        "count",
        "--config",
        "nonexistent.toml",
        "--input-file",
        "test_input.csv",
        "--output-file",
        "test_output.csv",
        "--barcode-column",
        "umi3_umi5_corrected",
        "--key-column",
        "protospacer",
    ]
    cli = Cli(args)
    with pytest.raises(ValueError, match="Configuration file not found"):
        cli.run()


def test_count_invalid_config_section():
    """Test that count command handles config file with invalid section"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file with wrong section name
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[wrong_section]
sep = ";"
"""
            )

        args = [
            "count",
            "--config",
            config_file,
            "--input-file",
            "test_input.csv",
            "--output-file",
            "test_output.csv",
            "--barcode-column",
            "umi3_umi5_corrected",
            "--key-column",
            "protospacer",
        ]
        cli = Cli(args)
        # Should still work, just using defaults for count-specific settings
        assert cli.args.sep == ","  # Default value

        
# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
