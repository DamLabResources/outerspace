"""Functional tests for the OUTERSPACE CLI workflow using iterative commands"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import os
import pytest
from outerspace.cli.main import Cli


def test_full_workflow(temp_workspace):
    """Test the full workflow from example.sh"""
    # Step 1: Run findseq commands for each sample
    pairs = [
        (
            "reads/409-4_S1_L002_R1_001.fastq.gz",
            "reads/409-4_S1_L002_R2_001.fastq.gz",
            "shuffle",
        ),
        (
            "reads/2-G1L9-M1_S9_L001_R1_001.fastq.gz",
            "reads/2-G1L9-M1_S9_L001_R2_001.fastq.gz",
            "M1-lib",
        ),
        (
            "reads/2-G1L9-M2_S12_L001_R1_001.fastq.gz",
            "reads/2-G1L9-M2_S12_L001_R2_001.fastq.gz",
            "M2-lib",
        ),
    ]

    for read1, read2, output_name in pairs:
        findseq_args = [
            "findseq",
            "-c", os.path.join(temp_workspace, "grnaquery.toml"),
            "-1",
            os.path.join(temp_workspace, read1),
            "-2",
            os.path.join(temp_workspace, read2),
            "-o",
            os.path.join(temp_workspace, f"findseq/{output_name}.csv"),
        ]
        cli = Cli(findseq_args)
        cli.run()

        # Verify findseq output
        assert os.path.exists(
            os.path.join(temp_workspace, f"findseq/{output_name}.csv")
        )

    # Step 2: Run collapse command
    collapse_args = [
        "collapse",
        "--input-dir",
        os.path.join(temp_workspace, "findseq"),
        "--output-dir",
        os.path.join(temp_workspace, "collapse"),
        "--columns",
        "UMI_5prime,UMI_3prime",
        "--mismatches",
        "2",
        "--method",
        "directional",
    ]
    cli = Cli(collapse_args)
    cli.run()

    # Verify collapse output
    for output_name in ["shuffle", "M1-lib", "M2-lib"]:
        assert os.path.exists(
            os.path.join(temp_workspace, f"collapse/{output_name}.csv")
        )

    # Step 3: Run count command
    count_args = [
        "count",
        "--input-dir",
        os.path.join(temp_workspace, "collapse"),
        "--output-dir",
        os.path.join(temp_workspace, "count"),
        "--barcode-column",
        "UMI_5prime_UMI_3prime_corrected",
        "--key-column",
        "protospacer",
    ]
    cli = Cli(count_args)
    cli.run()

    # Verify count output
    for output_name in ["shuffle", "M1-lib", "M2-lib"]:
        assert os.path.exists(os.path.join(temp_workspace, f"count/{output_name}.csv"))

    # Step 4: Run merge command
    merge_args = [
        "merge",
        os.path.join(temp_workspace, "count/shuffle.csv"),
        os.path.join(temp_workspace, "count/M1-lib.csv"),
        os.path.join(temp_workspace, "count/M2-lib.csv"),
        "--output-file",
        os.path.join(temp_workspace, "merged_counts.csv"),
        "--key-column",
        "protospacer",
        "--count-column",
        "UMI_5prime_UMI_3prime_corrected_count",
        "--sample-names",
        "shuffle",
        "M1-lib",
        "M2-lib",
        "--format",
        "wide",
    ]
    cli = Cli(merge_args)
    cli.run()

    # Verify merge output
    assert os.path.exists(os.path.join(temp_workspace, "merged_counts.csv"))

    # Step 5: Run stats command
    stats_args = [
        "stats",
        os.path.join(temp_workspace, "count/shuffle.csv"),
        os.path.join(temp_workspace, "count/M1-lib.csv"),
        os.path.join(temp_workspace, "count/M2-lib.csv"),
        "--key-column",
        "protospacer",
        "--count-column",
        "UMI_5prime_UMI_3prime_corrected_count",
    ]
    cli = Cli(stats_args)
    cli.run()


def test_workflow_with_allowed_list(temp_workspace):
    """Test the workflow with an allowed list for counting"""
    # Create allowed list file
    allowed_list_path = os.path.join(temp_workspace, "library_protospacers.txt")

    # Run the workflow steps
    # Step 1: findseq (same as before)
    pairs = [
        (
            "reads/409-4_S1_L002_R1_001.fastq.gz",
            "reads/409-4_S1_L002_R2_001.fastq.gz",
            "shuffle",
        ),
        (
            "reads/2-G1L9-M1_S9_L001_R1_001.fastq.gz",
            "reads/2-G1L9-M1_S9_L001_R2_001.fastq.gz",
            "M1-lib",
        ),
        (
            "reads/2-G1L9-M2_S12_L001_R1_001.fastq.gz",
            "reads/2-G1L9-M2_S12_L001_R2_001.fastq.gz",
            "M2-lib",
        ),
    ]

    for read1, read2, output_name in pairs:
        findseq_args = [
            "findseq",
            "-c", os.path.join(temp_workspace, "grnaquery.toml"),
            "-1",
            os.path.join(temp_workspace, read1),
            "-2",
            os.path.join(temp_workspace, read2),
            "-o",
            os.path.join(temp_workspace, f"findseq/{output_name}.csv"),
        ]
        cli = Cli(findseq_args)
        cli.run()

    # Step 2: collapse (same as before)
    collapse_args = [
        "collapse",
        "--input-dir",
        os.path.join(temp_workspace, "findseq"),
        "--output-dir",
        os.path.join(temp_workspace, "collapse"),
        "--columns",
        "UMI_5prime,UMI_3prime",
        "--mismatches",
        "2",
        "--method",
        "directional",
    ]
    cli = Cli(collapse_args)
    cli.run()

    # Step 3: count with allowed list
    count_args = [
        "count",
        "--input-dir",
        os.path.join(temp_workspace, "collapse"),
        "--output-dir",
        os.path.join(temp_workspace, "count"),
        "--barcode-column",
        "UMI_5prime_UMI_3prime_corrected",
        "--key-column",
        "protospacer",
        "--allowed-list",
        allowed_list_path,
    ]
    cli = Cli(count_args)
    cli.run()

    # Verify count output
    for output_name in ["shuffle", "M1-lib", "M2-lib"]:
        assert os.path.exists(os.path.join(temp_workspace, f"count/{output_name}.csv"))

    # Step 4: Run merge command with long format
    merge_args = [
        "merge",
        os.path.join(temp_workspace, "count/shuffle.csv"),
        os.path.join(temp_workspace, "count/M1-lib.csv"),
        os.path.join(temp_workspace, "count/M2-lib.csv"),
        "--output-file",
        os.path.join(temp_workspace, "merged_counts_long.csv"),
        "--key-column",
        "protospacer",
        "--count-column",
        "UMI_5prime_UMI_3prime_corrected_count",
        "--sample-names",
        "shuffle",
        "M1-lib",
        "M2-lib",
        "--format",
        "long",
    ]
    cli = Cli(merge_args)
    cli.run()

    # Verify merge output
    assert os.path.exists(os.path.join(temp_workspace, "merged_counts_long.csv"))

    # Step 5: Run stats command
    stats_args = [
        "stats",
        os.path.join(temp_workspace, "count/shuffle.csv"),
        os.path.join(temp_workspace, "count/M1-lib.csv"),
        os.path.join(temp_workspace, "count/M2-lib.csv"),
        "--key-column",
        "protospacer",
        "--count-column",
        "UMI_5prime_UMI_3prime_corrected_count",
    ]
    cli = Cli(stats_args)
    cli.run()


def test_full_workflow_with_key_rescue(temp_workspace):
    """Run the full iterative workflow but enable key-rescue during count.

    This ensures the commands interoperate without errors when key rescue is on.
    """
    # Step 1: Run findseq commands for each sample
    pairs = [
        (
            "reads/409-4_S1_L002_R1_001.fastq.gz",
            "reads/409-4_S1_L002_R2_001.fastq.gz",
            "shuffle",
        ),
        (
            "reads/2-G1L9-M1_S9_L001_R1_001.fastq.gz",
            "reads/2-G1L9-M1_S9_L001_R2_001.fastq.gz",
            "M1-lib",
        ),
        (
            "reads/2-G1L9-M2_S12_L001_R1_001.fastq.gz",
            "reads/2-G1L9-M2_S12_L001_R2_001.fastq.gz",
            "M2-lib",
        ),
    ]

    for read1, read2, output_name in pairs:
        findseq_args = [
            "findseq",
            "-c", os.path.join(temp_workspace, "grnaquery.toml"),
            "-1",
            os.path.join(temp_workspace, read1),
            "-2",
            os.path.join(temp_workspace, read2),
            "-o",
            os.path.join(temp_workspace, f"findseq/{output_name}.csv"),
        ]
        cli = Cli(findseq_args)
        cli.run()

    # Step 2: Run collapse command
    collapse_args = [
        "collapse",
        "--input-dir",
        os.path.join(temp_workspace, "findseq"),
        "--output-dir",
        os.path.join(temp_workspace, "collapse"),
        "--columns",
        "UMI_5prime,UMI_3prime",
        "--mismatches",
        "2",
        "--method",
        "directional",
    ]
    cli = Cli(collapse_args)
    cli.run()

    # Step 3: Run count command with allowed list and key rescue
    count_args = [
        "count",
        "--input-dir",
        os.path.join(temp_workspace, "collapse"),
        "--output-dir",
        os.path.join(temp_workspace, "count"),
        "--barcode-column",
        "UMI_5prime_UMI_3prime_corrected",
        "--key-column",
        "protospacer",
        "--allowed-list",
        os.path.join(temp_workspace, "library_protospacers.txt"),
        "--key-rescue",
        "--key-mismatch-penalty", "-1",
        "--key-gap-penalty", "-3",
        "--key-match-score", "1",
        "--key-min-score", "2",
    ]
    cli = Cli(count_args)
    cli.run()

    # Verify count output files exist
    for output_name in ["shuffle", "M1-lib", "M2-lib"]:
        assert os.path.exists(os.path.join(temp_workspace, f"count/{output_name}.csv"))

    # Step 4: Run merge command
    merge_args = [
        "merge",
        os.path.join(temp_workspace, "count/shuffle.csv"),
        os.path.join(temp_workspace, "count/M1-lib.csv"),
        os.path.join(temp_workspace, "count/M2-lib.csv"),
        "--output-file",
        os.path.join(temp_workspace, "merged_counts.csv"),
        "--key-column",
        "protospacer",
        "--count-column",
        "UMI_5prime_UMI_3prime_corrected_count",
        "--sample-names",
        "shuffle",
        "M1-lib",
        "M2-lib",
        "--format",
        "wide",
    ]
    cli = Cli(merge_args)
    cli.run()

    assert os.path.exists(os.path.join(temp_workspace, "merged_counts.csv"))

    # Step 5: Run stats command
    stats_args = [
        "stats",
        os.path.join(temp_workspace, "count/shuffle.csv"),
        os.path.join(temp_workspace, "count/M1-lib.csv"),
        os.path.join(temp_workspace, "count/M2-lib.csv"),
        "--key-column",
        "protospacer",
        "--count-column",
        "UMI_5prime_UMI_3prime_corrected_count",
    ]
    cli = Cli(stats_args)
    cli.run()

def test_single_file_workflow(temp_workspace):
    """Test the workflow using single file mode for collapse and count"""
    # Step 1: Run findseq for a single sample
    findseq_args = [
        "findseq",
        "-c", os.path.join(temp_workspace, "grnaquery.toml"),
        "-1",
        os.path.join(temp_workspace, "reads/409-4_S1_L002_R1_001.fastq.gz"),
        "-2",
        os.path.join(temp_workspace, "reads/409-4_S1_L002_R2_001.fastq.gz"),
        "-o",
        os.path.join(temp_workspace, "findseq/shuffle.csv"),
    ]
    cli = Cli(findseq_args)
    cli.run()

    # Verify findseq output
    assert os.path.exists(os.path.join(temp_workspace, "findseq/shuffle.csv"))

    # Step 2: Run collapse command on single file
    collapse_args = [
        "collapse",
        "--input-file",
        os.path.join(temp_workspace, "findseq/shuffle.csv"),
        "--output-file",
        os.path.join(temp_workspace, "collapse/shuffle.csv"),
        "--columns",
        "UMI_5prime,UMI_3prime",
        "--mismatches",
        "2",
        "--method",
        "directional",
    ]
    cli = Cli(collapse_args)
    cli.run()

    # Verify collapse output
    assert os.path.exists(os.path.join(temp_workspace, "collapse/shuffle.csv"))

    # Step 3: Run count command on single file
    count_args = [
        "count",
        "--input-file",
        os.path.join(temp_workspace, "collapse/shuffle.csv"),
        "--output-file",
        os.path.join(temp_workspace, "count/shuffle.csv"),
        "--barcode-column",
        "UMI_5prime_UMI_3prime_corrected",
        "--key-column",
        "protospacer",
    ]
    cli = Cli(count_args)
    cli.run()

    # Verify count output
    assert os.path.exists(os.path.join(temp_workspace, "count/shuffle.csv"))

    # Step 4: Run stats command on single file
    stats_args = [
        "stats",
        os.path.join(temp_workspace, "count/shuffle.csv"),
        "--key-column",
        "protospacer",
        "--count-column",
        "UMI_5prime_UMI_3prime_corrected_count",
    ]
    cli = Cli(stats_args)
    cli.run()


def test_workflow_with_config(temp_workspace):
    """Test the full workflow using config file for options"""
    # Step 1: Run findseq commands for each sample
    pairs = [
        (
            "reads/409-4_S1_L002_R1_001.fastq.gz",
            "reads/409-4_S1_L002_R2_001.fastq.gz",
            "shuffle",
        ),
        (
            "reads/2-G1L9-M1_S9_L001_R1_001.fastq.gz",
            "reads/2-G1L9-M1_S9_L001_R2_001.fastq.gz",
            "M1-lib",
        ),
        (
            "reads/2-G1L9-M2_S12_L001_R1_001.fastq.gz",
            "reads/2-G1L9-M2_S12_L001_R2_001.fastq.gz",
            "M2-lib",
        ),
    ]

    for read1, read2, output_name in pairs:
        findseq_args = [
            "findseq",
            "-c", os.path.join(temp_workspace, "grnaquery.toml"),
            "-1",
            os.path.join(temp_workspace, read1),
            "-2",
            os.path.join(temp_workspace, read2),
            "-o",
            os.path.join(temp_workspace, f"findseq/{output_name}.csv"),
        ]
        cli = Cli(findseq_args)
        cli.run()

        # Verify findseq output
        assert os.path.exists(
            os.path.join(temp_workspace, f"findseq/{output_name}.csv")
        )

    # Step 2: Run collapse command using config
    collapse_args = [
        "collapse",
        "--input-dir",
        os.path.join(temp_workspace, "findseq"),
        "--output-dir",
        os.path.join(temp_workspace, "collapse"),
        "--config",
        os.path.join(temp_workspace, "grnaquery.toml"),
    ]
    cli = Cli(collapse_args)
    cli.run()

    # Verify collapse output
    for output_name in ["shuffle", "M1-lib", "M2-lib"]:
        assert os.path.exists(
            os.path.join(temp_workspace, f"collapse/{output_name}.csv")
        )

    # Step 3: Run count command using config
    count_args = [
        "count",
        "--input-dir",
        os.path.join(temp_workspace, "collapse"),
        "--output-dir",
        os.path.join(temp_workspace, "count"),
        "--config",
        os.path.join(temp_workspace, "grnaquery.toml"),
    ]
    cli = Cli(count_args)
    cli.run()

    # Verify count output
    for output_name in ["shuffle", "M1-lib", "M2-lib"]:
        assert os.path.exists(os.path.join(temp_workspace, f"count/{output_name}.csv"))

    # Step 4: Run merge command using config
    merge_args = [
        "merge",
        os.path.join(temp_workspace, "count/shuffle.csv"),
        os.path.join(temp_workspace, "count/M1-lib.csv"),
        os.path.join(temp_workspace, "count/M2-lib.csv"),
        "--output-file",
        os.path.join(temp_workspace, "merged_counts.csv"),
        "--config",
        os.path.join(temp_workspace, "grnaquery.toml"),
        "--sample-names",
        "shuffle",
        "M1-lib",
        "M2-lib",
    ]
    cli = Cli(merge_args)
    cli.run()

    # Verify merge output
    assert os.path.exists(os.path.join(temp_workspace, "merged_counts.csv"))

    # Step 5: Run stats command using config
    stats_args = [
        "stats",
        os.path.join(temp_workspace, "count/shuffle.csv"),
        os.path.join(temp_workspace, "count/M1-lib.csv"),
        os.path.join(temp_workspace, "count/M2-lib.csv"),
        "--config",
        os.path.join(temp_workspace, "grnaquery.toml"),
    ]
    cli = Cli(stats_args)
    cli.run()

    
# New: key-rescue functional test using a synthetic collapse file
def test_count_key_rescue_synthetic(temp_workspace):
    """Count should rescue near-miss keys to allowed keys with --key-rescue"""
    # Prepare synthetic collapse directory and file
    collapse_dir = os.path.join(temp_workspace, "collapse_rescue")
    os.makedirs(collapse_dir, exist_ok=True)
    collapse_file = os.path.join(collapse_dir, "synthetic.csv")

    with open(collapse_file, "w") as f:
        f.write("protospacer,UMI_5prime_UMI_3prime_corrected\n")
        f.write("key1,bar1\n")
        f.write("key1,bar2\n")
        f.write("kay1,bar3\n")  # 1 mismatch from key1

    # Allowed list containing only key1
    allowed_list = os.path.join(temp_workspace, "allowed_keys.txt")
    with open(allowed_list, "w") as f:
        f.write("key1\n")

    # Run count with key rescue enabled
    output_file = os.path.join(temp_workspace, "count_rescue.csv")
    args = [
        "count",
        "--input-file",
        collapse_file,
        "--output-file",
        output_file,
        "--barcode-column",
        "UMI_5prime_UMI_3prime_corrected",
        "--key-column",
        "protospacer",
        "--allowed-list",
        allowed_list,
        "--key-rescue",
        "--key-mismatch-penalty",
        "-1",
        "--key-gap-penalty",
        "-3",
        "--key-match-score",
        "1",
        "--key-min-score",
        "2",
    ]
    cli = Cli(args)
    cli.run()

    # Verify rescued result
    assert os.path.exists(output_file)
    with open(output_file, "r") as f:
        content = f.read()
        assert "key1,3" in content  # kay1 rescued to key1
        assert "kay1" not in content

# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.