"""Tests for stepwise stats functionality with _from_step methods"""

import csv
import os
import pytest
import tempfile
from pathlib import Path

from outerspace.stats import (
    GiniCoefficient,
    ShannonDiversity,
    SimpsonDiversity,
    UMIRecoveryRate,
    UMIEfficiencyRate,
    ErrorRate,
)

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


@pytest.fixture
def simple_csv_file():
    """Create a simple CSV file for testing"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
        writer = csv.DictWriter(f, fieldnames=['key', 'barcode', 'count'])
        writer.writeheader()
        writer.writerow({'key': 'AAAA', 'barcode': 'B1', 'count': '10'})
        writer.writerow({'key': 'AAAA', 'barcode': 'B2', 'count': '5'})
        writer.writerow({'key': 'TTTT', 'barcode': 'B1', 'count': '8'})
        writer.writerow({'key': 'CCCC', 'barcode': 'B3', 'count': '12'})
        writer.writerow({'key': 'GGGG', 'barcode': 'B1', 'count': '3'})
        filepath = f.name
    
    yield filepath
    os.unlink(filepath)


@pytest.fixture
def error_rate_csv_file():
    """Create a CSV file with original and corrected columns"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
        writer = csv.DictWriter(f, fieldnames=['original', 'corrected'])
        writer.writeheader()
        writer.writerow({'original': 'AAAA', 'corrected': 'AAAA'})  # No error
        writer.writerow({'original': 'ATTT', 'corrected': 'AAAA'})  # 3 errors
        writer.writerow({'original': 'AAAA', 'corrected': 'AAAA'})  # No error
        writer.writerow({'original': 'AAAT', 'corrected': 'AAAA'})  # 1 error
        writer.writerow({'original': 'TTTT', 'corrected': 'TTTT'})  # No error
        filepath = f.name
    
    yield filepath
    os.unlink(filepath)


@pytest.fixture
def allowed_list_file():
    """Create an allowed list file"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
        f.write("AAAA\n")
        f.write("TTTT\n")
        f.write("CCCC\n")
        f.write("GGGG\n")
        f.write("MISSING1\n")
        f.write("MISSING2\n")
        filepath = f.name
    
    yield filepath
    os.unlink(filepath)


def test_gini_from_step(simple_csv_file):
    """Test GiniCoefficient._from_step"""
    result = GiniCoefficient._from_step(
        simple_csv_file,
        key_column="key",
        barcode_column="count"
    )
    
    assert result is not None
    assert 0 <= result <= 1


def test_gini_from_step_with_allowed_list(simple_csv_file, allowed_list_file):
    """Test GiniCoefficient._from_step with allowed list"""
    result = GiniCoefficient._from_step(
        simple_csv_file,
        key_column="key",
        barcode_column="count",
        allowed_list=allowed_list_file
    )
    
    assert result is not None
    assert 0 <= result <= 1


def test_shannon_from_step(simple_csv_file):
    """Test ShannonDiversity._from_step"""
    result = ShannonDiversity._from_step(
        simple_csv_file,
        key_column="key",
        barcode_column="count"
    )
    
    assert result is not None
    assert result >= 0


def test_shannon_from_step_with_base(simple_csv_file):
    """Test ShannonDiversity._from_step with custom base"""
    result_base2 = ShannonDiversity._from_step(
        simple_csv_file,
        key_column="key",
        barcode_column="count",
        base=2.0
    )
    
    result_base_e = ShannonDiversity._from_step(
        simple_csv_file,
        key_column="key",
        barcode_column="count",
        base=2.718281828
    )
    
    assert result_base2 is not None
    assert result_base_e is not None
    # Results should be different with different bases
    assert abs(result_base2 - result_base_e) > 0.01


def test_simpson_from_step(simple_csv_file):
    """Test SimpsonDiversity._from_step"""
    result = SimpsonDiversity._from_step(
        simple_csv_file,
        key_column="key",
        barcode_column="count"
    )
    
    assert result is not None
    assert 0 <= result <= 1


def test_recovery_rate_from_step(simple_csv_file, allowed_list_file):
    """Test UMIRecoveryRate._from_step with allowed list"""
    result = UMIRecoveryRate._from_step(
        simple_csv_file,
        key_column="key",
        barcode_column="count",
        allowed_list=allowed_list_file
    )
    
    assert result is not None
    assert 0 <= result <= 1
    # We have 4 out of 6 keys in allowed list
    assert result == pytest.approx(4/6, rel=0.01)


def test_efficiency_rate_from_step(simple_csv_file, allowed_list_file):
    """Test UMIEfficiencyRate._from_step with allowed list"""
    result = UMIEfficiencyRate._from_step(
        simple_csv_file,
        key_column="key",
        barcode_column="count",
        allowed_list=allowed_list_file
    )
    
    assert result is not None
    assert 0 <= result <= 1


def test_error_rate_from_step(error_rate_csv_file):
    """Test ErrorRate._from_step with alignment scoring"""
    result = ErrorRate._from_step(
        error_rate_csv_file,
        original_column="original",
        corrected_column="corrected"
    )
    
    assert result is not None
    # Error rate should be between 0 and 1
    assert 0 <= result <= 1
    # With alignment scoring, error rate should reflect mismatches
    # The test file has some errors, so it should be > 0
    assert result > 0


def test_error_rate_missing_columns(simple_csv_file):
    """Test ErrorRate._from_step with missing columns"""
    with pytest.raises(ValueError, match="not found"):
        ErrorRate._from_step(
            simple_csv_file,
            original_column="nonexistent",
            corrected_column="corrected"
        )


def test_missing_required_params(simple_csv_file):
    """Test that missing required parameters raise errors"""
    with pytest.raises(ValueError, match="key_column is required"):
        GiniCoefficient._from_step(simple_csv_file)
    
    with pytest.raises(ValueError, match="key_column is required"):
        ShannonDiversity._from_step(simple_csv_file)


def test_recovery_rate_missing_allowed_list(simple_csv_file):
    """Test that UMIRecoveryRate requires allowed_list"""
    with pytest.raises(ValueError, match="allowed_list is required"):
        UMIRecoveryRate._from_step(
            simple_csv_file,
            key_column="key"
        )


def test_efficiency_rate_missing_allowed_list(simple_csv_file):
    """Test that UMIEfficiencyRate requires allowed_list"""
    with pytest.raises(ValueError, match="allowed_list is required"):
        UMIEfficiencyRate._from_step(
            simple_csv_file,
            key_column="key"
        )


def test_error_rate_missing_params(simple_csv_file):
    """Test that ErrorRate requires both columns"""
    with pytest.raises(ValueError, match="Both original_column and corrected_column are required"):
        ErrorRate._from_step(
            simple_csv_file,
            original_column="original"
        )


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.

