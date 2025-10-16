"""Tests for Shannon diversity calculations"""

import pytest
import numpy as np
from outerspace.stats import ShannonDiversity

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


def test_shannon_perfect_equality():
    """Test Shannon diversity calculation for perfectly equal distribution"""
    counts = [10, 10, 10, 10, 10]
    shannon = ShannonDiversity.calculate_shannon(counts)
    # For 5 equally abundant items, H = log2(5) ≈ 2.32
    assert abs(shannon - 2.32) < 0.01


def test_shannon_perfect_inequality():
    """Test Shannon diversity calculation for perfectly unequal distribution"""
    counts = [1e10, 1, 1, 1, 1]
    shannon = ShannonDiversity.calculate_shannon(counts)
    assert abs(shannon) < 0.0001  # Should be very close to 0


def test_shannon_moderate_inequality():
    """Test Shannon diversity calculation for moderately unequal distribution"""
    counts = [50, 30, 20, 10, 5]
    shannon = ShannonDiversity.calculate_shannon(counts)
    assert 0 < shannon < 2.32  # Should be between 0 and log2(5)


def test_shannon_different_bases():
    """Test Shannon diversity calculation with different logarithm bases"""
    counts = [10, 10, 10, 10, 10]
    shannon_base2 = ShannonDiversity.calculate_shannon(counts, base=2.0)
    shannon_base10 = ShannonDiversity.calculate_shannon(counts, base=10.0)
    assert shannon_base2 != shannon_base10
    assert abs(shannon_base2 / shannon_base10 - np.log2(10)) < 0.01


def test_shannon_empty_input():
    """Test Shannon diversity calculation with empty input"""
    shannon = ShannonDiversity.calculate_shannon([])
    assert shannon is None


def test_shannon_zero_counts():
    """Test Shannon diversity calculation with all zero counts"""
    counts = [0, 0, 0, 0, 0]
    shannon = ShannonDiversity.calculate_shannon(counts)
    assert shannon is None


def test_shannon_mixed_zero_counts():
    """Test Shannon diversity calculation with some zero counts"""
    counts = [10, 0, 10, 0, 10]
    shannon = ShannonDiversity.calculate_shannon(counts)
    # For 3 equally abundant items, H = log2(3) ≈ 1.58
    assert abs(shannon - 1.58) < 0.01


def test_shannon_end_to_end(equal_umi):
    """Test end-to-end Shannon diversity calculation"""
    result = ShannonDiversity.calculate(equal_umi)
    assert result is not None
    assert result > 0

    
# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.