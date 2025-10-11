"""Tests for the collapse command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
import os
import tempfile
from outerspace.cli.main import Cli


def test_collapse_initialization():
    """Test that collapse command initializes correctly"""
    args = [
        "collapse",
        "--input-dir",
        "test_input",
        "--output-dir",
        "test_output",
        "--columns",
        "umi3,umi5",
        "--mismatches",
        "2",
        "--method",
        "directional",
    ]
    cli = Cli(args)
    assert cli.args.input_dir == "test_input"
    assert cli.args.output_dir == "test_output"
    assert cli.args.columns == "umi3,umi5"
    assert cli.args.mismatches == 2
    assert cli.args.method == "directional"
    assert cli.args.sep == ","  # default value


def test_collapse_single_file_initialization():
    """Test that collapse command initializes correctly in single file mode"""
    args = [
        "collapse",
        "--input-file",
        "test_input.csv",
        "--output-file",
        "test_output.csv",
        "--columns",
        "umi3,umi5",
        "--mismatches",
        "2",
        "--method",
        "directional",
    ]
    cli = Cli(args)
    assert cli.args.input_file == "test_input.csv"
    assert cli.args.output_file == "test_output.csv"
    assert cli.args.columns == "umi3,umi5"
    assert cli.args.mismatches == 2
    assert cli.args.method == "directional"
    assert cli.args.sep == ","  # default value


def test_collapse_missing_input():
    """Test that collapse command handles missing input"""
    args = [
        "collapse",
        "--output-dir",
        "test_output",
        "--columns",
        "umi3,umi5",
        "--mismatches",
        "2",
        "--method",
        "directional",
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2


def test_collapse_missing_output():
    """Test that collapse command handles missing output"""
    args = [
        "collapse",
        "--input-dir",
        "test_input",
        "--columns",
        "umi3,umi5",
        "--mismatches",
        "2",
        "--method",
        "directional",
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2


def test_collapse_missing_columns():
    """Test that collapse command handles missing columns"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[collapse]
columns = "UMI_5prime,UMI_3prime"
mismatches = 2
method = "directional"
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("header1,header2\nvalue1,value2\n")

        args = [
            "collapse",
            "--input-file",
            input_file,
            "--output-file",
            os.path.join(temp_dir, "output.csv"),
            "--mismatches",
            "2",
            "--method",
            "directional",
        ]
        with pytest.raises(
            ValueError, match="Please provide either --columns or --config"
        ):
            cli = Cli(args)
            cli.run()


def test_collapse_invalid_method():
    """Test that collapse command handles invalid clustering method"""
    args = [
        "collapse",
        "--input-dir",
        "test_input",
        "--output-dir",
        "test_output",
        "--columns",
        "umi3,umi5",
        "--mismatches",
        "2",
        "--method",
        "invalid_method",
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2


def test_collapse_nonexistent_input_dir():
    """Test that collapse command handles nonexistent input directory"""
    args = [
        "collapse",
        "--input-dir",
        "nonexistent_dir",
        "--output-dir",
        "test_output",
        "--columns",
        "umi3,umi5",
        "--mismatches",
        "2",
        "--method",
        "directional",
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()


def test_collapse_empty_input_dir():
    """Test that collapse command handles empty input directory"""
    with tempfile.TemporaryDirectory() as temp_dir:
        args = [
            "collapse",
            "--input-dir",
            temp_dir,
            "--output-dir",
            "test_output",
            "--columns",
            "umi3,umi5",
            "--mismatches",
            "2",
            "--method",
            "directional",
        ]
        cli = Cli(args)
        with pytest.raises(ValueError):
            cli.run()


def test_collapse_parse_columns():
    """Test that column parsing works correctly"""
    args = [
        "collapse",
        "--input-dir",
        "test_input",
        "--output-dir",
        "test_output",
        "--columns",
        "umi3,umi5",
        "--mismatches",
        "2",
        "--method",
        "directional",
    ]
    cli = Cli(args)
    columns = cli.command._parse_columns(cli.args.columns)
    assert columns == ["umi3", "umi5"]

    # Test with spaces
    cli.args.columns = "umi3, umi5"
    columns = cli.command._parse_columns(cli.args.columns)
    assert columns == ["umi3", "umi5"]

    # Test with single column
    cli.args.columns = "umi3"
    columns = cli.command._parse_columns(cli.args.columns)
    assert columns == ["umi3"]


def test_collapse_single_file_missing_output_file():
    """Test that collapse command handles missing output file for single file mode"""
    args = [
        "collapse",
        "--input-file",
        "test_input.csv",
        "--columns",
        "umi3,umi5",
        "--mismatches",
        "2",
        "--method",
        "directional",
    ]
    with pytest.raises(SystemExit) as excinfo:
        cli = Cli(args)
    assert excinfo.value.code == 2


def test_collapse_single_file_nonexistent_input():
    """Test that collapse command handles nonexistent input file in single file mode"""
    args = [
        "collapse",
        "--input-file",
        "nonexistent.csv",
        "--output-file",
        "test_output.csv",
        "--columns",
        "umi3,umi5",
        "--mismatches",
        "2",
        "--method",
        "directional",
    ]
    cli = Cli(args)
    with pytest.raises(ValueError):
        cli.run()


def test_collapse_single_file_basic():
    """Test basic single file collapse functionality"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("umi3,umi5\n")
            f.write("AAA,CCC\n")
            f.write("AAT,CCC\n")  # One mismatch from AAA
            f.write("GGG,TTT\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "collapse",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--columns",
            "umi3,umi5",
            "--mismatches",
            "1",
            "--method",
            "directional",
        ]
        cli = Cli(args)
        cli.run()

        # Verify output file exists and has correct content
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            content = f.read()
            assert "umi3,umi5,umi3_umi5_corrected" in content
            # The first two rows should have the same corrected barcode
            assert "AAA,CCC" in content
            assert "AAT,CCC" in content
            assert "GGG,TTT" in content


def test_collapse_with_config_file():
    """Test that collapse command loads and uses config file correctly"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[collapse]
method = "cluster"
mismatches = 3
sep = ";"
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("umi3;umi5\n")
            f.write("AAA;CCC\n")
            f.write("AAT;CCC\n")
            f.write("GGG;TTT\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "collapse",
            "--config",
            config_file,
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--columns",
            "umi3,umi5",
        ]
        cli = Cli(args)
        cli.run()

        # Verify config values were applied
        assert cli.args.method == "cluster"
        assert cli.args.mismatches == 3
        assert cli.args.sep == ";"

        # Verify output file exists and has correct content
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            content = f.read()
            assert "umi3;umi5;umi3_umi5_corrected" in content


def test_collapse_config_override():
    """Test that command line arguments override config file settings"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[collapse]
method = "cluster"
mismatches = 3
sep = ";"
"""
            )

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("umi3;umi5\n")  # Using semicolon separator to match config
            f.write("AAA;CCC\n")
            f.write("AAT;CCC\n")
            f.write("GGG;TTT\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "collapse",
            "--config",
            config_file,
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--columns",
            "umi3,umi5",
            "--method",
            "adjacency",  # Override config
            "--mismatches",
            "1",  # Override config
        ]
        cli = Cli(args)
        cli.run()

        # Verify command line args took precedence
        assert cli.args.method == "adjacency"
        assert cli.args.mismatches == 1
        assert cli.args.sep == ";"  # This one should come from config

        # Verify output file exists and has correct content
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            content = f.read()
            assert "umi3;umi5;umi3_umi5_corrected" in content


def test_collapse_nonexistent_config():
    """Test that collapse command handles nonexistent config file"""
    args = [
        "collapse",
        "--config",
        "nonexistent.toml",
        "--input-file",
        "test_input.csv",
        "--output-file",
        "test_output.csv",
        "--columns",
        "umi3,umi5",
    ]
    cli = Cli(args)
    with pytest.raises(ValueError, match="Configuration file not found"):
        cli.run()


def test_collapse_invalid_config_section():
    """Test that collapse command handles config file with invalid section"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file with wrong section name
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write(
                """[wrong_section]
method = "cluster"
"""
            )

        args = [
            "collapse",
            "--config",
            config_file,
            "--input-file",
            "test_input.csv",
            "--output-file",
            "test_output.csv",
            "--columns",
            "umi3,umi5",
        ]
        cli = Cli(args)
        # Should still work, just using defaults for collapse-specific settings
        assert cli.args.method == "directional"  # Default value


def test_collapse_with_nearest_method():
    """Test collapse with method=nearest (performs key rescue)"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # First, create a file with UMI-corrected data
        umi_corrected_file = os.path.join(temp_dir, "umi_corrected.csv")
        with open(umi_corrected_file, "w") as f:
            f.write("protospacer,umi3_umi5_corrected\n")
            f.write("key1,AAACCC\n")
            f.write("key1,AAACCC\n")
            f.write("kayX,GGGTT\n")  # Near-miss for key1
            f.write("other,TTTAAA\n")

        # Create allowed list file
        allowed_list = os.path.join(temp_dir, "allowed.txt")
        with open(allowed_list, "w") as f:
            f.write("key1\n")
            f.write("other\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "proto_corrected.csv")

        # Run collapse with method=nearest
        args = [
            "collapse",
            "--input-file",
            umi_corrected_file,
            "--output-file",
            output_file,
            "--columns",
            "protospacer",
            "--allowed-list",
            allowed_list,
            "--rescue-exhaustive",
            "--mismatch-penalty",
            "-1",
            "--gap-penalty",
            "-3",
            "--match-score",
            "1",
            "--min-score",
            "0",
            "--method",
            "nearest",
        ]
        cli = Cli(args)
        cli.run()

        # Verify output file exists and has corrected column
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            content = f.read()
            assert "protospacer_corrected" in content
            # kayX should be rescued to key1
            lines = content.strip().split("\n")
            # Find the row with kayX
            for line in lines[1:]:  # Skip header
                if "kayX" in line:
                    # Check that the corrected column has key1
                    assert "key1" in line


def test_collapse_with_allowed_method():
    """Test collapse with method=allowed (exact matching only, no rescue)"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file (simulating UMI-corrected data)
        input_file = os.path.join(temp_dir, "umi_corrected.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,umi3_umi5_corrected\n")
            f.write("key1,AAACCC\n")
            f.write("kayX,GGGTT\n")  # Near key1 but should be filtered (not rescued)
            f.write("key2,TTTAAA\n")

        # Create allowed list file
        allowed_list = os.path.join(temp_dir, "allowed.txt")
        with open(allowed_list, "w") as f:
            f.write("key1\n")
            f.write("key2\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "proto_corrected.csv")

        args = [
            "collapse",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--columns",
            "protospacer",
            "--allowed-list",
            allowed_list,
            "--method",
            "allowed",
        ]
        cli = Cli(args)
        cli.run()

        # Verify output file exists and has filtered keys
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            content = f.read()
            assert "protospacer_corrected" in content
            lines = content.strip().split("\n")
            # Find the row with kayX
            for line in lines[1:]:  # Skip header
                if "kayX" in line:
                    # Check that the corrected column is empty
                    fields = line.split(",")
                    # The last field should be the protospacer_corrected (empty)
                    assert fields[-1] == ""


def test_collapse_nearest_method_requires_allowed_list():
    """Test that method=nearest requires allowed list"""
    with tempfile.TemporaryDirectory() as temp_dir:
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,umi3_umi5_corrected\n")
            f.write("key1,AAACCC\n")

        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "collapse",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--columns",
            "protospacer",
            "--method",
            "nearest",
            # No allowed list
        ]
        cli = Cli(args)
        with pytest.raises(ValueError, match="Method 'nearest' requires --allowed-list"):
            cli.run()


def test_collapse_allowed_method_requires_allowed_list():
    """Test that method=allowed requires allowed list"""
    with tempfile.TemporaryDirectory() as temp_dir:
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,umi3_umi5_corrected\n")
            f.write("key1,AAACCC\n")

        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "collapse",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--columns",
            "protospacer",
            "--method",
            "allowed",
            # No allowed list
        ]
        cli = Cli(args)
        with pytest.raises(ValueError, match="Method 'allowed' requires --allowed-list"):
            cli.run()


def test_collapse_with_steps():
    """Test collapse with iterative steps from config"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create input directory with test files
        input_dir = os.path.join(temp_dir, "findseq")
        os.makedirs(input_dir)
        
        input_file = os.path.join(input_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,UMI_5prime,UMI_3prime\n")
            f.write("key1,AAAA,CCCC\n")
            f.write("key1,AAAT,CCCC\n")  # 1 mismatch in UMI
            f.write("key2,GGGG,TTTT\n")  # Near-miss for key1 (1 mismatch)
            f.write("other,TTTT,AAAA\n")

        # Create allowed list file
        allowed_list = os.path.join(temp_dir, "allowed.txt")
        with open(allowed_list, "w") as f:
            f.write("key1\n")
            f.write("other\n")

        # Create config file with steps
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write("""
[[collapse.steps]]
name = "umi_correction"
columns = "UMI_5prime,UMI_3prime"
method = "directional"
mismatches = 1

[[collapse.steps]]
name = "protospacer_correction"
columns = "protospacer"
method = "nearest"
allowed_list = "allowed.txt"
mismatch_penalty = -1
gap_penalty = -3
match_score = 1
min_score = 0
""")

        # Create output directory
        output_dir = os.path.join(temp_dir, "output")

        args = [
            "collapse",
            "--config",
            config_file,
            "--input-dir",
            input_dir,
            "--output-dir",
            output_dir,
        ]
        cli = Cli(args)
        cli.run()

        # Verify output file exists
        output_file = os.path.join(output_dir, "test.csv")
        assert os.path.exists(output_file)
        
        with open(output_file, "r") as f:
            content = f.read()
            # Should have both corrected columns
            assert "UMI_5prime_UMI_3prime_corrected" in content
            assert "protospacer_corrected" in content
            # key2 should be rescued to key1 (1 mismatch, score = 2*1 + 1*(-1) = 1)
            lines = content.strip().split("\n")
            rescued = False
            for line in lines[1:]:  # Skip header
                if "key2," in line:
                    # Check that the corrected column has key1
                    parts = line.split(",")
                    # Last column should be protospacer_corrected
                    if len(parts) > 4:
                        assert parts[-1] == "key1", f"Expected key2 to be rescued to key1, but got {parts[-1]}"
                        rescued = True
            assert rescued, "key2 should have been found and rescued to key1"


def test_collapse_with_steps_single_file():
    """Test collapse with iterative steps from config in single file mode"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create single input file (not in a directory)
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,UMI_5prime,UMI_3prime\n")
            f.write("key1,AAAA,CCCC\n")
            f.write("key1,AAAT,CCCC\n")  # 1 mismatch in UMI
            f.write("key2,GGGG,TTTT\n")  # Near-miss for key1 (1 mismatch)
            f.write("other,TTTT,AAAA\n")

        # Create allowed list file
        allowed_list = os.path.join(temp_dir, "allowed.txt")
        with open(allowed_list, "w") as f:
            f.write("key1\n")
            f.write("other\n")

        # Create config file with steps
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write("""
[[collapse.steps]]
name = "umi_correction"
columns = "UMI_5prime,UMI_3prime"
method = "directional"
mismatches = 1

[[collapse.steps]]
name = "protospacer_correction"
columns = "protospacer"
method = "nearest"
allowed_list = "allowed.txt"
mismatch_penalty = -1
gap_penalty = -3
match_score = 1
min_score = 0
""")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "collapse",
            "--config",
            config_file,
            "--input-file",
            input_file,
            "--output-file",
            output_file,
        ]
        cli = Cli(args)
        cli.run()

        # Verify output file exists
        assert os.path.exists(output_file)
        
        with open(output_file, "r") as f:
            content = f.read()
            # Should have both corrected columns
            assert "UMI_5prime_UMI_3prime_corrected" in content
            assert "protospacer_corrected" in content
            # key2 should be rescued to key1
            lines = content.strip().split("\n")
            rescued = False
            for line in lines[1:]:  # Skip header
                if "key2," in line:
                    # Check that the corrected column has key1
                    parts = line.split(",")
                    # Last column should be protospacer_corrected
                    if len(parts) > 4:
                        assert parts[-1] == "key1", f"Expected key2 to be rescued to key1, but got {parts[-1]}"
                        rescued = True
            assert rescued, "key2 should have been found and rescued to key1"


# ============================================================================
# Multithreading Tests
# ============================================================================

def test_collapse_with_threads_parameter():
    """Test collapse command with --threads parameter"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,umi\n")
            f.write("key1,AAA\n")
            f.write("key1,AAA\n")
            f.write("kayX,BBB\n")  # Near-miss for key1
            f.write("other,CCC\n")

        # Create allowed list file
        allowed_list = os.path.join(temp_dir, "allowed.txt")
        with open(allowed_list, "w") as f:
            f.write("key1\n")
            f.write("other\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "collapse",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--columns",
            "protospacer",
            "--allowed-list",
            allowed_list,
            "--method",
            "nearest",
            "--min-score",
            "0",
            "--threads",
            "2",
        ]
        cli = Cli(args)
        
        # Verify threads parameter is set
        assert cli.args.threads == 2
        
        cli.run()

        # Verify output file exists and has corrected column
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            content = f.read()
            assert "protospacer_corrected" in content


def test_collapse_threads_default_value():
    """Test that threads parameter defaults to 1"""
    args = [
        "collapse",
        "--input-file",
        "test.csv",
        "--output-file",
        "output.csv",
        "--columns",
        "umi",
    ]
    cli = Cli(args)
    assert cli.args.threads == 1


def test_collapse_with_threads_consistency():
    """Test that results are consistent with different thread counts"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file with multiple rows
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,umi\n")
            for i in range(20):
                f.write(f"key{i % 3},UMI{i:02d}\n")
            # Add some near-misses
            f.write("kay1,UMI99\n")  # Near key1
            f.write("kay2,UMI98\n")  # Near key2

        # Create allowed list file
        allowed_list = os.path.join(temp_dir, "allowed.txt")
        with open(allowed_list, "w") as f:
            f.write("key0\n")
            f.write("key1\n")
            f.write("key2\n")

        # Run with threads=1
        output_file_1 = os.path.join(temp_dir, "output1.csv")
        args_1 = [
            "collapse",
            "--input-file",
            input_file,
            "--output-file",
            output_file_1,
            "--columns",
            "protospacer",
            "--allowed-list",
            allowed_list,
            "--method",
            "nearest",
            "--min-score",
            "0",
            "--threads",
            "1",
        ]
        cli_1 = Cli(args_1)
        cli_1.run()

        # Run with threads=4
        output_file_4 = os.path.join(temp_dir, "output4.csv")
        args_4 = [
            "collapse",
            "--input-file",
            input_file,
            "--output-file",
            output_file_4,
            "--columns",
            "protospacer",
            "--allowed-list",
            allowed_list,
            "--method",
            "nearest",
            "--min-score",
            "0",
            "--threads",
            "4",
        ]
        cli_4 = Cli(args_4)
        cli_4.run()

        # Read both outputs and compare
        with open(output_file_1, "r") as f:
            content_1 = f.read()
        with open(output_file_4, "r") as f:
            content_4 = f.read()

        # Results should be identical
        assert content_1 == content_4


def test_collapse_with_threads_in_config():
    """Test that threads parameter can be specified in config file"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file with threads
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write("""[collapse]
threads = 3
method = "nearest"
min_score = 0
""")

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,umi\n")
            f.write("key1,AAA\n")
            f.write("other,BBB\n")

        # Create allowed list file
        allowed_list = os.path.join(temp_dir, "allowed.txt")
        with open(allowed_list, "w") as f:
            f.write("key1\n")
            f.write("other\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "collapse",
            "--config",
            config_file,
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--columns",
            "protospacer",
            "--allowed-list",
            allowed_list,
        ]
        cli = Cli(args)
        cli.run()
        
        # Verify config was loaded
        assert cli.args.threads == 3

        # Verify output file exists
        assert os.path.exists(output_file)


def test_collapse_threads_cli_overrides_config():
    """Test that CLI --threads parameter overrides config file"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file with threads=2
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write("""[collapse]
threads = 2
method = "nearest"
min_score = 0
""")

        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,umi\n")
            f.write("key1,AAA\n")

        # Create allowed list file
        allowed_list = os.path.join(temp_dir, "allowed.txt")
        with open(allowed_list, "w") as f:
            f.write("key1\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "collapse",
            "--config",
            config_file,
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--columns",
            "protospacer",
            "--allowed-list",
            allowed_list,
            "--threads",
            "5",  # Override config
        ]
        cli = Cli(args)
        
        # Verify CLI arg took precedence
        assert cli.args.threads == 5


def test_collapse_with_threads_in_steps():
    """Test collapse with threads parameter in steps mode"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("protospacer,umi\n")
            for i in range(10):
                f.write(f"key{i % 2},UMI{i:02d}\n")
            f.write("kay1,UMI99\n")  # Near key1

        # Create allowed list file
        allowed_list = os.path.join(temp_dir, "allowed.txt")
        with open(allowed_list, "w") as f:
            f.write("key0\n")
            f.write("key1\n")

        # Create config file with steps and threads
        config_file = os.path.join(temp_dir, "config.toml")
        with open(config_file, "w") as f:
            f.write("""
[[collapse.steps]]
name = "protospacer_correction"
columns = "protospacer"
method = "nearest"
allowed_list = "allowed.txt"
threads = 3
min_score = 0
""")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        args = [
            "collapse",
            "--config",
            config_file,
            "--input-file",
            input_file,
            "--output-file",
            output_file,
        ]
        cli = Cli(args)
        cli.run()

        # Verify output file exists
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            content = f.read()
            assert "protospacer_corrected" in content


def test_collapse_threads_with_directory_mode():
    """Test collapse with threads in directory mode"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create input directory with multiple files
        input_dir = os.path.join(temp_dir, "input")
        os.makedirs(input_dir)
        
        for file_num in range(3):
            input_file = os.path.join(input_dir, f"test{file_num}.csv")
            with open(input_file, "w") as f:
                f.write("protospacer,umi\n")
                f.write(f"key{file_num},AAA\n")
                f.write(f"kay{file_num},BBB\n")  # Near-miss

        # Create allowed list file
        allowed_list = os.path.join(temp_dir, "allowed.txt")
        with open(allowed_list, "w") as f:
            f.write("key0\n")
            f.write("key1\n")
            f.write("key2\n")

        # Create output directory
        output_dir = os.path.join(temp_dir, "output")

        args = [
            "collapse",
            "--input-dir",
            input_dir,
            "--output-dir",
            output_dir,
            "--columns",
            "protospacer",
            "--allowed-list",
            allowed_list,
            "--method",
            "nearest",
            "--min-score",
            "0",
            "--threads",
            "2",
        ]
        cli = Cli(args)
        cli.run()

        # Verify output files exist
        for file_num in range(3):
            output_file = os.path.join(output_dir, f"test{file_num}.csv")
            assert os.path.exists(output_file)


def test_collapse_threads_only_affects_nearest_method():
    """Test that threads parameter is only used with method=nearest"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test input file
        input_file = os.path.join(temp_dir, "test.csv")
        with open(input_file, "w") as f:
            f.write("umi3,umi5\n")
            f.write("AAA,CCC\n")
            f.write("AAT,CCC\n")  # One mismatch from AAA
            f.write("GGG,TTT\n")

        # Create output file path
        output_file = os.path.join(temp_dir, "output.csv")

        # Test with directional method (threads shouldn't affect it)
        args = [
            "collapse",
            "--input-file",
            input_file,
            "--output-file",
            output_file,
            "--columns",
            "umi3,umi5",
            "--mismatches",
            "1",
            "--method",
            "directional",
            "--threads",
            "4",  # This should be ignored for directional method
        ]
        cli = Cli(args)
        cli.run()

        # Verify output file exists
        assert os.path.exists(output_file)
        with open(output_file, "r") as f:
            content = f.read()
            assert "umi3_umi5_corrected" in content


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
