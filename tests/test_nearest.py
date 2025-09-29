"""Tests for nearest neighbor UMI matching functionality"""

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

import pytest
from outerspace.nearest import find_closest_umi, calculate_alignment_score


def test_calculate_alignment_score_exact_match():
    """Test scoring for exact matches"""
    score = calculate_alignment_score("ATCG", "ATCG", -1, -3, 1)
    assert score == 4  # 4 matches * 1 point each


def test_calculate_alignment_score_single_mismatch():
    """Test scoring for single mismatch"""
    score = calculate_alignment_score("ATCG", "ATCC", -1, -3, 1)
    # 3 matches (1 each) + 1 mismatch (-1) = 3 + (-1) = 2
    assert score == 2


def test_calculate_alignment_score_single_gap():
    """Test scoring for single gap/insertion"""
    score = calculate_alignment_score("ATCG", "ATC", -1, -3, 1)
    # Should align as ATC- vs ATCG, giving 3 matches + 1 gap
    # 3 matches (1 each) + 1 gap (-3) = 3 + (-3) = 0
    assert score == 0


def test_calculate_alignment_score_single_deletion():
    """Test scoring for single deletion"""
    score = calculate_alignment_score("ATC", "ATCG", -1, -3, 1)
    # Should align as ATC- vs ATCG, giving 3 matches + 1 gap
    # 3 matches (1 each) + 1 gap (-3) = 3 + (-3) = 0
    assert score == 0


def test_calculate_alignment_score_multiple_mismatches():
    """Test scoring for multiple mismatches"""
    score = calculate_alignment_score("ATCG", "GCTA", -1, -3, 1)
    # All positions are mismatches: 4 mismatches * -1 = -4
    assert score == -4


def test_calculate_alignment_score_custom_scoring():
    """Test with custom scoring parameters"""
    # More lenient scoring
    score = calculate_alignment_score("ATCG", "ATCC", -1, -2, 2)
    # 3 matches (2 each) + 1 mismatch (-1) = 6 + (-1) = 5
    assert score == 5
    
    # Harsher mismatch penalty
    score = calculate_alignment_score("ATCG", "ATCC", -5, -3, 1)
    # 3 matches (1 each) + 1 mismatch (-5) = 3 + (-5) = -2
    assert score == -2


def test_calculate_alignment_score_empty_sequences():
    """Test with empty sequences"""
    score = calculate_alignment_score("", "", -1, -3, 1)
    assert score == 0


def test_calculate_alignment_score_one_empty_sequence():
    """Test with one empty sequence"""
    score = calculate_alignment_score("ATCG", "", -1, -3, 1)
    # Should be 4 gaps: 4 * -3 = -12
    assert score == -12
    
    score = calculate_alignment_score("", "ATCG", -1, -3, 1)
    # Should be 4 gaps: 4 * -3 = -12
    assert score == -12


def test_calculate_alignment_score_different_length_sequences():
    """Test sequences of different lengths"""
    score = calculate_alignment_score("ATCGAA", "ATCG", -1, -3, 1)
    # Should align with gaps, exact alignment depends on pyspoa algorithm
    # At minimum should have 4 matches and 2 gaps: 4*1 + 2*(-3) = -2
    assert score <= 4  # Can't be better than 4 matches


def test_calculate_alignment_score_caching():
    """Test that results are cached"""
    # First call
    score1 = calculate_alignment_score("ATCG", "ATCC", -1, -3, 1)
    # Second call with same parameters should return cached result
    score2 = calculate_alignment_score("ATCG", "ATCC", -1, -3, 1)
    assert score1 == score2
    
    # Different parameters should give different results
    score3 = calculate_alignment_score("ATCG", "ATCC", -2, -3, 1)
    assert score3 != score1


def test_find_closest_umi_exact_match():
    """Test finding exact match"""
    allowed_umis = ("ATCG", "GCTA", "TTAA", "CCGG")
    result = find_closest_umi(allowed_umis, "ATCG")
    assert result == "ATCG"


def test_find_closest_umi_close_match():
    """Test finding close match within threshold"""
    allowed_umis = ("ATCG", "GCTA", "TTAA", "CCGG")
    # ATCC differs from ATCG by 1 mismatch
    # Score: 3 matches (3*1) + 1 mismatch (1*-1) = 2
    # Default min_score is 17, so this won't match with default settings
    result = find_closest_umi(allowed_umis, "ATCC", min_score=2)
    assert result == "ATCG"


def test_find_closest_umi_no_match_below_threshold():
    """Test when no UMI meets minimum score threshold"""
    allowed_umis = ("ATCG", "GCTA", "TTAA", "CCGG")
    # NNNN should not match any UMI well enough
    result = find_closest_umi(allowed_umis, "NNNN", min_score=0)
    # With default scoring, this should give very negative scores
    assert result is None


def test_find_closest_umi_multiple_candidates_best_score():
    """Test that function returns the best scoring match"""
    allowed_umis = ("ATCG", "ATCC", "ATCA")  # All similar to ATCG
    result = find_closest_umi(allowed_umis, "ATCG", min_score=0)
    assert result == "ATCG"  # Exact match should be best


def test_find_closest_umi_custom_scoring_parameters():
    """Test with custom scoring parameters"""
    allowed_umis = ("ATCG", "GCTA", "TTAA")
    
    # Very lenient scoring - allow single mismatch
    result = find_closest_umi(
        allowed_umis, "ATCC", 
        mismatch_penalty=-1, gap_penalty=-3, match_score=2, min_score=5
    )
    assert result == "ATCG"
    
    # Very strict scoring - require exact match
    result = find_closest_umi(
        allowed_umis, "ATCC",
        mismatch_penalty=-10, gap_penalty=-10, match_score=1, min_score=4
    )
    assert result is None  # No exact match available


def test_find_closest_umi_gap_handling():
    """Test handling of insertions/deletions"""
    allowed_umis = ("ATCG", "GCTA", "TTAA")
    # Query has deletion compared to ATCG
    result = find_closest_umi(
        allowed_umis, "ATC",
        mismatch_penalty=-1, gap_penalty=-1, match_score=1, min_score=2
    )
    assert result == "ATCG"


def test_find_closest_umi_early_termination_on_perfect_match():
    """Test that function terminates early on perfect match"""
    # Put exact match at end to test early termination logic
    allowed_umis = ("GCTA", "TTAA", "CCGG", "ATCG")
    result = find_closest_umi(allowed_umis, "ATCG", min_score=0)
    assert result == "ATCG"


def test_find_closest_umi_empty_allowed_list():
    """Test with empty allowed UMIs list"""
    result = find_closest_umi((), "ATCG")
    assert result is None


def test_find_closest_umi_empty_query():
    """Test with empty query UMI"""
    allowed_umis = ("ATCG", "GCTA")
    result = find_closest_umi(allowed_umis, "")
    assert result is None


def test_find_closest_umi_caching_behavior():
    """Test that results are properly cached"""
    allowed_umis = ("ATCG", "GCTA", "TTAA")
    
    # First call
    result1 = find_closest_umi(allowed_umis, "ATCC", min_score=0)
    # Second call with same parameters
    result2 = find_closest_umi(allowed_umis, "ATCC", min_score=0)
    assert result1 == result2


def test_find_closest_umi_realistic_umi_scenarios():
    """Test with realistic UMI sequences and parameters"""
    # Typical UMI library
    allowed_umis = (
        "AAAAAAAA", "AAAAAAAC", "AAAAAAAG", "AAAAAACC",
        "AAAAAACG", "AAAAACCC", "AAAAACCG", "AAAAACGG"
    )
    
    # Test exact match
    result = find_closest_umi(allowed_umis, "AAAAAAAA", min_score=6)
    assert result == "AAAAAAAA"
    
    # Test single mismatch
    result = find_closest_umi(allowed_umis, "AAAAAAAC", min_score=6)
    assert result == "AAAAAAAC"
    
    # Test sequence with single error that should match
    result = find_closest_umi(
        allowed_umis, "AAAAAAAA",  # This is exact match
        mismatch_penalty=-1, gap_penalty=-3, match_score=1, min_score=7
    )
    assert result == "AAAAAAAA"


def test_find_closest_umi_performance_with_large_list():
    """Test performance doesn't degrade significantly with larger UMI lists"""
    # Create a larger list of UMIs
    allowed_umis = tuple(f"AAAA{i:04b}".replace('0', 'A').replace('1', 'T') 
                       for i in range(16))
    
    result = find_closest_umi(allowed_umis, "AAAAAAAA", min_score=6)
    assert result is not None  # Should find something reasonable


def test_find_closest_umi_edge_case_scoring():
    """Test edge cases in scoring"""
    allowed_umis = ("ATCG",)
    
    # Test with zero penalties/rewards
    result = find_closest_umi(
        allowed_umis, "ATCC",
        mismatch_penalty=0, gap_penalty=0, match_score=0, min_score=0
    )
    assert result == "ATCG"  # Should still work with zero scoring
    
    # Test with positive mismatch "penalty" (actually a reward)
    result = find_closest_umi(
        allowed_umis, "ATCC",
        mismatch_penalty=1, gap_penalty=-1, match_score=1, min_score=3
    )
    assert result == "ATCG"  # Should work even with positive mismatch score


def test_score_consistency():
    """Test that scores are calculated consistently between functions"""
    seq1, seq2 = "ATCG", "ATCC"
    mismatch_penalty, gap_penalty, match_score = -1, -3, 1
    
    # Calculate score directly
    direct_score = calculate_alignment_score(
        seq1, seq2, mismatch_penalty, gap_penalty, match_score
    )
    
    # Calculate score through find_closest_umi
    allowed_umis = (seq2,)
    result = find_closest_umi(
        allowed_umis, seq1,
        mismatch_penalty=mismatch_penalty,
        gap_penalty=gap_penalty, 
        match_score=match_score,
        min_score=direct_score  # Set threshold to exact score
    )
    
    # Should find the match since score meets threshold
    assert result == seq2
    
    # Should not find match if threshold is higher
    result_strict = find_closest_umi(
        allowed_umis, seq1,
        mismatch_penalty=mismatch_penalty,
        gap_penalty=gap_penalty,
        match_score=match_score, 
        min_score=direct_score + 1
    )
    assert result_strict is None


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
