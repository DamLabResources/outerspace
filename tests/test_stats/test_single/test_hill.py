"""Tests for Hill number calculations"""

import pytest
import numpy as np
from outerspace.stats import HillNumber

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


def test_hill_q0_richness():
    """Test Hill number with q=0 (species richness)"""
    counts = {
        b"AAAAAA": 10,
        b"TTTTTT": 20,
        b"CCCCCC": 30,
        b"GGGGGG": 5,
    }
    result = HillNumber.calculate_hill(list(counts.values()), q=0)
    # q=0 should return number of species
    assert result == 4


def test_hill_q1_shannon():
    """Test Hill number with q=1 (exponential Shannon)"""
    counts = {
        b"AAAAAA": 10,
        b"TTTTTT": 10,
        b"CCCCCC": 10,
    }
    result = HillNumber.calculate_hill(list(counts.values()), q=1)
    # q=1 with equal counts should equal number of species
    assert abs(result - 3.0) < 0.01


def test_hill_q2_simpson():
    """Test Hill number with q=2 (inverse Simpson)"""
    counts = {
        b"AAAAAA": 10,
        b"TTTTTT": 20,
        b"CCCCCC": 30,
    }
    result = HillNumber.calculate_hill(list(counts.values()), q=2)
    # q=2 is inverse Simpson concentration
    # proportions: 1/6, 1/3, 1/2
    # sum(p^2) = 1/36 + 1/9 + 1/4 = 1/36 + 4/36 + 9/36 = 14/36
    # D_2 = 1 / (14/36) = 36/14 = 2.571...
    expected = 1 / ((10/60)**2 + (20/60)**2 + (30/60)**2)
    assert abs(result - expected) < 0.01


def test_hill_general_q():
    """Test Hill number with arbitrary q value"""
    counts = {
        b"AAAAAA": 10,
        b"TTTTTT": 20,
        b"CCCCCC": 30,
    }
    result = HillNumber.calculate_hill(list(counts.values()), q=1.5)
    assert result is not None
    assert result > 0


def test_hill_empty_counts():
    """Test Hill number with empty counts"""
    counts = {}
    result = HillNumber.calculate_hill(list(counts.values()), q=1)
    assert result is None


def test_hill_zero_counts():
    """Test Hill number with zero counts"""
    counts = {
        b"AAAAAA": 0,
        b"TTTTTT": 0,
    }
    result = HillNumber.calculate_hill(list(counts.values()), q=1)
    assert result is None


def test_hill_with_zeros():
    """Test Hill number with some zero counts (should be filtered out)"""
    counts = {
        b"AAAAAA": 10,
        b"TTTTTT": 0,
        b"CCCCCC": 20,
    }
    result = HillNumber.calculate_hill(list(counts.values()), q=0)
    # Should only count non-zero species
    assert result == 2


def test_hill_end_to_end(equal_umi):
    """Test end-to-end Hill number calculation"""
    result = HillNumber.calculate(equal_umi, q=1)
    assert result is not None
    assert result > 0


def test_hill_with_allowed_list(partial_umi):
    """Test Hill number calculation with allowed list"""
    allowed_list = ["AAAAAA", "TTTTTT", "CCCCCC"]
    result = HillNumber.calculate(partial_umi, q=1, allowed_list=allowed_list)
    assert result is not None
    assert result > 0


def test_hill_different_q_values():
    """Test that Hill numbers are ordered: q=0 >= q=1 >= q=2"""
    counts = {
        b"AAAAAA": 10,
        b"TTTTTT": 20,
        b"CCCCCC": 30,
        b"GGGGGG": 5,
    }
    
    q0 = HillNumber.calculate_hill(list(counts.values()), q=0)
    q1 = HillNumber.calculate_hill(list(counts.values()), q=1)
    q2 = HillNumber.calculate_hill(list(counts.values()), q=2)
    
    # Hill numbers should decrease as q increases (Renyi diversity)
    assert q0 >= q1
    assert q1 >= q2


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.

