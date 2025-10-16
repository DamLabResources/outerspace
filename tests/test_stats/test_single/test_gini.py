"""Tests for Gini coefficient calculation"""

import pytest
from outerspace.stats import GiniCoefficient

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


def test_gini_perfect_equality(equal_umi):
    """Test Gini coefficient calculation for perfectly equal distribution"""
    result = GiniCoefficient.calculate(equal_umi)
    gini = result
    assert abs(gini) < 0.0001


def test_gini_perfect_inequality(unequal_umi):
    """Test Gini coefficient calculation for perfectly unequal distribution"""
    result = GiniCoefficient.calculate(unequal_umi)
    gini = result
    assert 0.8 < gini < 1.0


def test_gini_moderate_inequality(moderate_umi):
    """Test Gini coefficient calculation for moderately unequal distribution"""
    result = GiniCoefficient.calculate(moderate_umi)
    gini = result
    assert 0 <= gini < 1
    assert 0.2 < gini < 0.3


def test_gini_empty_input(empty_umi):
    """Test Gini coefficient calculation with empty input"""
    result = GiniCoefficient.calculate(empty_umi)
    assert result is None


def test_gini_zero_counts():
    """Test Gini coefficient calculation with all zero counts"""
    from outerspace.umi import UMI

    umi = UMI(mismatches=0)
    sequences = ["AAAAAA", "TTTTTT", "CCCCCC"]
    for seq in sequences:
        umi.consume(seq, 0)  # Add with zero count
    umi.create_mapping()

    result = GiniCoefficient.calculate(umi)
    assert result is None

    
# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.