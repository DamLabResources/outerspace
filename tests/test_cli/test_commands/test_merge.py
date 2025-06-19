"""Tests for the merge command"""

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
import os
import tempfile
import csv
import pandas as pd
from outerspace.cli.main import Cli

def test_merge_initialization():
    """Test that merge command initializes correctly"""
    args = [
        'merge',
        'test1.csv',
        'test2.csv',
        '--output-file', 'test_output.csv',
        '--key-column', 'umi',
        '--count-column', 'count'
    ]
    cli = Cli(args)
    assert cli.args.files == ['test1.csv', 'test2.csv']
    assert cli.args.output_file == 'test_output.csv'
    assert cli.args.key_column == 'umi'
    assert cli.args.count_column == 'count'
    assert cli.args.format == 'wide'  # default format

def test_merge_missing_files():
    """Test that merge command handles missing input files"""
    args = [
        'merge',
        'nonexistent1.csv',
        'nonexistent2.csv',
        '--output-file', 'test_output.csv',
        '--key-column', 'umi'
    ]
    cli = Cli(args)
    with pytest.raises(ValueError, match="Input file not found"):
        cli.run()

def test_merge_missing_key_column():
    """Test that merge command requires key column"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create config file
        config_file = os.path.join(temp_dir, 'config.toml')
        with open(config_file, 'w') as f:
            f.write("""[merge]
column = "UMI_5prime_UMI_3prime_corrected_count"
key_column = "protospacer"
""")

        # Create test input files
        input_files = []
        for i in range(2):
            input_file = os.path.join(temp_dir, f'test{i}.csv')
            with open(input_file, 'w') as f:
                f.write('header1,header2\nvalue1,value2\n')
            input_files.append(input_file)

        args = [
            'merge',
            *input_files,
            '--output-file', os.path.join(temp_dir, 'output.csv')
        ]
        with pytest.raises(ValueError, match="Please provide either --key-column or --config"):
            cli = Cli(args)
            cli.run()

def test_merge_sample_names_mismatch():
    """Test that merge command handles sample names mismatch"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test files
        file1 = os.path.join(temp_dir, 'test1.csv')
        file2 = os.path.join(temp_dir, 'test2.csv')
        
        # Write test data
        with open(file1, 'w') as f:
            f.write('umi,count\nACTG,10\nGCTA,5\n')
        with open(file2, 'w') as f:
            f.write('umi,count\nACTG,8\nGCTA,12\n')
        
        args = [
            'merge',
            file1,
            file2,
            '--output-file', os.path.join(temp_dir, 'output.csv'),
            '--key-column', 'umi',
            '--sample-names', 'sample1'  # Only one name for two files
        ]
        cli = Cli(args)
        with pytest.raises(ValueError, match="Number of sample names must match number of files"):
            cli.run()

def test_merge_with_example_data():
    """Test merge command with real example data"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test files
        file1 = os.path.join(temp_dir, 'test1.csv')
        file2 = os.path.join(temp_dir, 'test2.csv')
        output_file = os.path.join(temp_dir, 'merged.csv')
        
        # Write test data
        with open(file1, 'w') as f:
            f.write('umi,count\nACTG,10\nGCTA,5\n')
        with open(file2, 'w') as f:
            f.write('umi,count\nACTG,8\nGCTA,12\n')
        
        args = [
            'merge',
            file1,
            file2,
            '--output-file', output_file,
            '--key-column', 'umi',
            '--count-column', 'count',
            '--sample-names', 'sample1', 'sample2'
        ]
        cli = Cli(args)
        cli.run()
        
        # Verify output file exists and has correct format
        assert os.path.exists(output_file)
        
        # Read output and verify contents
        df = pd.read_csv(output_file, index_col='umi')
        assert set(df.columns) == {'sample1', 'sample2'}
        assert df.loc['ACTG', 'sample1'] == 10
        assert df.loc['ACTG', 'sample2'] == 8
        assert df.loc['GCTA', 'sample1'] == 5
        assert df.loc['GCTA', 'sample2'] == 12

def test_merge_long_format():
    """Test merge command with long format output"""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create test files
        file1 = os.path.join(temp_dir, 'test1.csv')
        file2 = os.path.join(temp_dir, 'test2.csv')
        output_file = os.path.join(temp_dir, 'merged.csv')
        
        # Write test data
        with open(file1, 'w') as f:
            f.write('umi,count\nACTG,10\nGCTA,5\n')
        with open(file2, 'w') as f:
            f.write('umi,count\nACTG,8\nGCTA,12\n')
        
        args = [
            'merge',
            file1,
            file2,
            '--output-file', output_file,
            '--key-column', 'umi',
            '--count-column', 'count',
            '--sample-names', 'sample1', 'sample2',
            '--format', 'long'
        ]
        cli = Cli(args)
        cli.run()
        
        # Verify output file exists and has correct format
        assert os.path.exists(output_file)
        
        # Read output and verify contents
        df = pd.read_csv(output_file)
        assert set(df.columns) == {'sample', 'umi', 'count'}
        assert len(df) == 4  # 2 UMIs * 2 samples
        
        # Verify specific values
        sample1_actg = df[(df['sample'] == 'sample1') & (df['umi'] == 'ACTG')]['count'].iloc[0]
        sample2_actg = df[(df['sample'] == 'sample2') & (df['umi'] == 'ACTG')]['count'].iloc[0]
        sample1_gcta = df[(df['sample'] == 'sample1') & (df['umi'] == 'GCTA')]['count'].iloc[0]
        sample2_gcta = df[(df['sample'] == 'sample2') & (df['umi'] == 'GCTA')]['count'].iloc[0]
        
        assert sample1_actg == 10
        assert sample2_actg == 8
        assert sample1_gcta == 5
        assert sample2_gcta == 12 