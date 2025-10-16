"""Tests for Simpson diversity calculations"""

import pytest
import numpy as np
from outerspace.stats import SimpsonDiversity

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


def test_simpson_perfect_equality():
    """Test Simpson diversity calculation for perfectly equal distribution"""
    counts = [10, 10, 10, 10, 10]
    simpson = SimpsonDiversity.calculate_simpson(counts)
    # For 5 equally abundant items, D = 1 - (1/5) = 0.8
    assert abs(simpson - 0.8) < 0.01


def test_simpson_perfect_inequality():
    """Test Simpson diversity calculation for perfectly unequal distribution"""
    counts = [100, 0, 0, 0, 0]
    simpson = SimpsonDiversity.calculate_simpson(counts)
    assert abs(simpson) < 0.0001  # Should be very close to 0


def test_simpson_moderate_inequality():
    """Test Simpson diversity calculation for moderately unequal distribution"""
    counts = [50, 30, 20, 10, 5]
    simpson = SimpsonDiversity.calculate_simpson(counts)
    assert 0 < simpson < 0.8  # Should be between 0 and 0.8


def test_simpson_empty_input():
    """Test Simpson diversity calculation with empty input"""
    simpson = SimpsonDiversity.calculate_simpson([])
    assert simpson is None


def test_simpson_zero_counts():
    """Test Simpson diversity calculation with all zero counts"""
    counts = [0, 0, 0, 0, 0]
    simpson = SimpsonDiversity.calculate_simpson(counts)
    assert simpson is None


def test_simpson_mixed_zero_counts():
    """Test Simpson diversity calculation with some zero counts"""
    counts = [10, 0, 10, 0, 10]
    simpson = SimpsonDiversity.calculate_simpson(counts)
    # For 3 equally abundant items, D = 1 - (1/3) ≈ 0.67
    assert abs(simpson - 0.67) < 0.01


def test_simpson_end_to_end(equal_umi):
    """Test end-to-end Simpson diversity calculation"""
    result = SimpsonDiversity.calculate(equal_umi)
    assert result is not None
    assert 0 < result < 1

    
# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.