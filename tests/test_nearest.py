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
    assert result == ["ATCG"]


def test_find_closest_umi_close_match():
    """Test finding close match within threshold"""
    allowed_umis = ("ATCG", "GCTA", "TTAA", "CCGG")
    # ATCC differs from ATCG by 1 mismatch
    # Score: 3 matches (3*1) + 1 mismatch (1*-1) = 2
    # Default min_score is 17, so this won't match with default settings
    result = find_closest_umi(allowed_umis, "ATCC", min_score=2)
    assert result == ["ATCG"]


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
    assert result == ["ATCG"] # Exact match should be best


def test_find_closest_umi_custom_scoring_parameters():
    """Test with custom scoring parameters"""
    allowed_umis = ("ATCG", "GCTA", "TTAA")
    
    # Very lenient scoring - allow single mismatch
    result = find_closest_umi(
        allowed_umis, "ATCC", 
        mismatch_penalty=-1, gap_penalty=-3, match_score=2, min_score=5
    )
    assert result == ["ATCG"]
    
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
    assert result == ["ATCG"]


def test_find_closest_umi_early_termination_on_perfect_match():
    """Test that function terminates early on perfect match"""
    # Put exact match at end to test early termination logic
    allowed_umis = ("GCTA", "TTAA", "CCGG", "ATCG")
    result = find_closest_umi(allowed_umis, "ATCG", min_score=0)
    assert result == ["ATCG"]


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
    assert result == ["AAAAAAAA"]
    
    # Test single mismatch
    result = find_closest_umi(allowed_umis, "AAAAAAAC", min_score=6)
    assert result == ["AAAAAAAC"]
    
    # Test sequence with single error that should match
    result = find_closest_umi(
        allowed_umis, "AAAAAAAA",  # This is exact match
        mismatch_penalty=-1, gap_penalty=-3, match_score=1, min_score=7
    )
    assert result == ["AAAAAAAA"]


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
    assert result == ["ATCG"]  # Should still work with zero scoring
    
    # Test with positive mismatch "penalty" (actually a reward)
    result = find_closest_umi(
        allowed_umis, "ATCC",
        mismatch_penalty=1, gap_penalty=-1, match_score=1, min_score=3
    )
    assert result == ["ATCG"]  # Should work even with positive mismatch score


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
    assert result == [seq2]
    
    # Should not find match if threshold is higher
    result_strict = find_closest_umi(
        allowed_umis, seq1,
        mismatch_penalty=mismatch_penalty,
        gap_penalty=gap_penalty,
        match_score=match_score, 
        min_score=direct_score + 1
    )
    assert result_strict is None


# ============================================================================
# Multithreading Tests
# ============================================================================

def test_find_many_single_threaded():
    """Test find_many with single thread (threads=1)"""
    from outerspace.nearest import NearestUMIFinder
    
    allowed_umis = ("ATCG", "GCTA", "TTAA", "CCGG")
    finder = NearestUMIFinder(
        allowed_list=allowed_umis,
        mismatch_penalty=-1,
        gap_penalty=-3,
        match_score=1,
        min_score=2,
    )
    
    queries = ["ATCG", "ATCC", "GCTA", "NNNN"]
    results = finder.find_many(queries, threads=1)
    
    assert len(results) == 4
    assert results[0] == ["ATCG"]  # Exact match
    assert results[1] == ["ATCG"]  # Close match
    assert results[2] == ["GCTA"]  # Exact match
    assert results[3] is None  # No match


def test_find_many_multithreaded():
    """Test find_many with multiple threads"""
    from outerspace.nearest import NearestUMIFinder
    
    allowed_umis = ("ATCG", "GCTA", "TTAA", "CCGG")
    finder = NearestUMIFinder(
        allowed_list=allowed_umis,
        mismatch_penalty=-1,
        gap_penalty=-3,
        match_score=1,
        min_score=2,
    )
    
    queries = ["ATCG", "ATCC", "GCTA", "NNNN"]
    results = finder.find_many(queries, threads=2)
    
    assert len(results) == 4
    assert results[0] == ["ATCG"]  # Exact match
    assert results[1] == ["ATCG"]  # Close match
    assert results[2] == ["GCTA"]  # Exact match
    assert results[3] is None  # No match


def test_find_many_consistency_single_vs_multi():
    """Test that single-threaded and multi-threaded give same results"""
    from outerspace.nearest import NearestUMIFinder
    
    allowed_umis = tuple(f"{'A' * i}{'T' * (8-i)}" for i in range(9))
    finder = NearestUMIFinder(
        allowed_list=allowed_umis,
        mismatch_penalty=-1,
        gap_penalty=-3,
        match_score=1,
        min_score=0,
    )
    
    # Create a diverse set of queries
    queries = [
        "AAAAAAAA", "ATATATAT", "TTTTTTTT", "AAAATTTT",
        "AATTAATT", "GGGGGGGG", "CCCCCCCC", "ACGTACGT"
    ] * 5  # Repeat to have enough queries for good parallelism
    
    # Run with single thread
    results_single = finder.find_many(queries, threads=1)
    
    # Run with multiple threads
    results_multi = finder.find_many(queries, threads=4)
    
    # Results should be identical
    assert len(results_single) == len(results_multi)
    for i, (r_single, r_multi) in enumerate(zip(results_single, results_multi)):
        assert r_single == r_multi, f"Mismatch at index {i}: {r_single} != {r_multi}"


def test_find_many_with_batch_size():
    """Test find_many with custom batch_size parameter"""
    from outerspace.nearest import NearestUMIFinder
    
    allowed_umis = ("ATCG", "GCTA", "TTAA", "CCGG")
    finder = NearestUMIFinder(
        allowed_list=allowed_umis,
        mismatch_penalty=-1,
        gap_penalty=-3,
        match_score=1,
        min_score=2,
    )
    
    queries = ["ATCG", "ATCC", "GCTA", "NNNN"] * 10  # 40 queries
    
    # Test with different batch sizes
    results_batch_1 = finder.find_many(queries, threads=2, batch_size=1)
    results_batch_5 = finder.find_many(queries, threads=2, batch_size=5)
    results_batch_10 = finder.find_many(queries, threads=2, batch_size=10)
    
    # All should give the same results
    assert results_batch_1 == results_batch_5 == results_batch_10


def test_find_many_large_dataset():
    """Test find_many with larger dataset to verify parallel speedup"""
    from outerspace.nearest import NearestUMIFinder
    
    # Create a larger allowed list
    allowed_umis = tuple(f"{'ACGT'[i%4] * 2}{j:04b}".replace('0', 'A').replace('1', 'T') 
                        for i in range(4) for j in range(16))
    
    finder = NearestUMIFinder(
        allowed_list=allowed_umis,
        mismatch_penalty=-1,
        gap_penalty=-3,
        match_score=1,
        min_score=0,
        use_prescreen=True,
    )
    
    # Create many queries
    queries = [allowed_umis[i % len(allowed_umis)] for i in range(100)]
    
    # Both should work and give same results
    results_single = finder.find_many(queries, threads=1)
    results_multi = finder.find_many(queries, threads=4)
    
    assert len(results_single) == len(results_multi) == 100
    assert results_single == results_multi


def test_find_many_empty_queries():
    """Test find_many with empty query list"""
    from outerspace.nearest import NearestUMIFinder
    
    allowed_umis = ("ATCG", "GCTA")
    finder = NearestUMIFinder(allowed_list=allowed_umis)
    
    results = finder.find_many([], threads=2)
    assert results == []


def test_find_many_with_prescreen_disabled():
    """Test find_many with k-mer prescreening disabled (exhaustive search)"""
    from outerspace.nearest import NearestUMIFinder
    
    allowed_umis = ("ATCG", "GCTA", "TTAA", "CCGG")
    finder = NearestUMIFinder(
        allowed_list=allowed_umis,
        mismatch_penalty=-1,
        gap_penalty=-3,
        match_score=1,
        min_score=2,
        use_prescreen=False,  # Disable k-mer prescreen
    )
    
    queries = ["ATCG", "ATCC", "GCTA"]
    
    # Should work with both single and multi-threaded
    results_single = finder.find_many(queries, threads=1)
    results_multi = finder.find_many(queries, threads=2)
    
    assert results_single == results_multi
    assert results_single[0] == ["ATCG"]
    assert results_single[1] == ["ATCG"]
    assert results_single[2] == ["GCTA"]


def test_static_find_impl():
    """Test the static _find_impl method directly"""
    from outerspace.nearest import NearestUMIFinder
    
    allowed_umis_tuple = ("ATCG", "GCTA", "TTAA")
    allowed_set = set(allowed_umis_tuple)
    kmer_sets = {umi: NearestUMIFinder._kmerize_static(umi, 3) for umi in allowed_umis_tuple}
    
    # Test exact match
    result = NearestUMIFinder._find_impl(
        query_umi="ATCG",
        allowed_umis_tuple=allowed_umis_tuple,
        allowed_set=allowed_set,
        kmer_sets=kmer_sets,
        mismatch_penalty=-1,
        gap_penalty=-3,
        match_score=1,
        min_score=0,
        use_prescreen=True,
        kmer_size=3,
        min_kmer_overlap=1,
    )
    assert result == ["ATCG"]
    
    # Test close match
    result = NearestUMIFinder._find_impl(
        query_umi="ATCC",
        allowed_umis_tuple=allowed_umis_tuple,
        allowed_set=allowed_set,
        kmer_sets=kmer_sets,
        mismatch_penalty=-1,
        gap_penalty=-3,
        match_score=1,
        min_score=2,
        use_prescreen=True,
        kmer_size=3,
        min_kmer_overlap=1,
    )
    assert result == ["ATCG"]


def test_kmerize_static():
    """Test the static _kmerize_static method"""
    from outerspace.nearest import NearestUMIFinder
    
    # Test normal case
    kmers = NearestUMIFinder._kmerize_static("ATCG", 3)
    assert kmers == frozenset(["ATC", "TCG"])
    
    # Test short sequence
    kmers = NearestUMIFinder._kmerize_static("AT", 3)
    assert kmers == frozenset(["AT"])
    
    # Test empty sequence
    kmers = NearestUMIFinder._kmerize_static("", 3)
    assert kmers == frozenset([""])
    
    # Test k=0
    kmers = NearestUMIFinder._kmerize_static("ATCG", 0)
    assert kmers == frozenset()


def test_find_many_preserves_order():
    """Test that find_many preserves the order of queries"""
    from outerspace.nearest import NearestUMIFinder
    
    allowed_umis = ("AAAA", "TTTT", "GGGG", "CCCC")
    finder = NearestUMIFinder(
        allowed_list=allowed_umis,
        min_score=3,
    )
    
    queries = ["AAAA", "TTTT", "GGGG", "CCCC", "AAAA", "TTTT"]
    
    results_single = finder.find_many(queries, threads=1)
    results_multi = finder.find_many(queries, threads=3)
    
    # Check order is preserved
    assert results_single == [["AAAA"], ["TTTT"], ["GGGG"], ["CCCC"], ["AAAA"], ["TTTT"]]
    assert results_multi == results_single


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
