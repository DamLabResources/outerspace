"""Tests for single UMI library statistics"""

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
from outerspace.umi import UMI, UmiCollection
from outerspace.stats.single import (
    UMIStats,
    GiniCoefficient,
    ShannonDiversity,
    SimpsonDiversity,
    HillNumber,
    UMIRecoveryRate,
    UMIEfficiencyRate,
    ErrorRate,
)

# List of all single UMI stat classes to test
SINGLE_STAT_CLASSES = [
    GiniCoefficient,
    ShannonDiversity,
    SimpsonDiversity,
    HillNumber,
    UMIRecoveryRate,
    UMIEfficiencyRate,
]

# Classes that support allowed_list
ALLOWED_LIST_CLASSES = [
    UMIRecoveryRate,
    UMIEfficiencyRate,
]


def test_single_stat_interface(empty_umi):
    """Test that single stats accept UMI objects and return correct types"""
    for stat_class in SINGLE_STAT_CLASSES:
        # Test calculate method
        result = stat_class.calculate(empty_umi)
        assert isinstance(
            result, (float, dict, type(None))
        ), f"{stat_class.__name__}.calculate() should return float, dict, or None"


def test_single_stat_collection_interface(empty_umi):
    """Test that single stats accept UmiCollection objects and return correct shapes"""
    # Create a collection with one UMI
    collection = UmiCollection({"sample1": empty_umi})

    for stat_class in SINGLE_STAT_CLASSES:
        # Test calculate_collection method
        results = stat_class.calculate_collection(collection)
        assert isinstance(
            results, dict
        ), f"{stat_class.__name__}.calculate_collection() should return a dict"
        assert (
            "sample1" in results
        ), f"{stat_class.__name__}.calculate_collection() should include all samples"
        assert isinstance(
            results["sample1"], (float, dict, type(None))
        ), f"{stat_class.__name__}.calculate_collection() values should be float, dict, or None"


def test_single_stat_collection_multiple_samples(empty_umi, equal_umi):
    """Test that single stats handle multiple samples correctly"""
    # Create a collection with multiple UMIs
    collection = UmiCollection({"sample1": empty_umi, "sample2": equal_umi})

    for stat_class in SINGLE_STAT_CLASSES:
        # Test calculate_collection method
        results = stat_class.calculate_collection(collection)
        assert isinstance(
            results, dict
        ), f"{stat_class.__name__}.calculate_collection() should return a dict"
        assert set(results.keys()) == {
            "sample1",
            "sample2",
        }, f"{stat_class.__name__}.calculate_collection() should include all samples"
        for value in results.values():
            assert isinstance(
                value, (float, dict, type(None))
            ), f"{stat_class.__name__}.calculate_collection() values should be float, dict, or None"


def test_single_stat_with_allowed_list(empty_umi):
    """Test that single stats handle allowed_list parameter correctly"""
    allowed_list = ["AAAAAA", "TTTTTT"]

    for stat_class in ALLOWED_LIST_CLASSES:
        # Test calculate method with allowed_list
        result = stat_class.calculate(empty_umi, allowed_list=allowed_list)
        assert isinstance(
            result, (float, dict, type(None))
        ), f"{stat_class.__name__}.calculate() should handle allowed_list parameter"

        # Test calculate_collection method with allowed_list
        collection = UmiCollection({"sample1": empty_umi})
        results = stat_class.calculate_collection(collection, allowed_list=allowed_list)
        assert isinstance(
            results, dict
        ), f"{stat_class.__name__}.calculate_collection() should handle allowed_list parameter"


# Add tests for other single-library statistics here


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.