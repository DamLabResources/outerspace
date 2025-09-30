"""Nearest neighbor matching functionality for UMIs.

This module provides functions for finding the closest matching UMI from a list
of allowed UMIs using a configurable global alignment scoring system.
"""

import logging
from typing import Optional, Tuple
from functools import lru_cache

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


@lru_cache(maxsize=1000)
def find_closest_umi(
    allowed_umis: Tuple[str, ...], 
    query_umi: str, 
    mismatch_penalty: int = -1,
    gap_penalty: int = -3,
    match_score: int = 1,
    min_score: int = 0
) -> Optional[str]:
    """Find the closest UMI match above a minimum score threshold.
    
    This function computes a global alignment score between the query and each
    candidate UMI using a simple Needleman–Wunsch style scoring scheme and
    returns the best scoring UMI if it meets the minimum score.
    
    Parameters
    ----------
    allowed_umis : Tuple[str, ...]
        Tuple of allowed UMI sequences to search against. Must be a tuple for hashing
        in the LRU cache.
    query_umi : str
        Query UMI sequence to find a match for
    mismatch_penalty : int, default=-1
        Penalty score for mismatches (should be negative)
    gap_penalty : int, default=-3
        Penalty score for insertions and deletions (should be negative)
    match_score : int, default=1
        Score for matches (should be positive)
    min_score : int, default=17
        Minimum alignment score required for a match. If no UMI scores above
        this threshold, returns None.
    
    Returns
    -------
    Optional[str]
        The best scoring UMI from allowed_umis if found above min_score,
        None otherwise.
        
    Higher scores indicate better matches. The scoring system allows for fine-tuned
    control over how mismatches and indels are penalized relative to matches.
    """
    if not allowed_umis or not query_umi:
        return None
    
    best_match = None
    best_score = min_score - 1  # Start below minimum threshold
    
    # Perfect score for an exact match of equal length
    max_possible_score = len(query_umi) * match_score

    # Check each allowed UMI for the best match
    for candidate_str in allowed_umis:
        
        # If exact match, short-circuit
        if candidate_str == query_umi:
            if max_possible_score >= min_score:
                return candidate_str
            # If an exact match does not meet the threshold, we still compute
            # other candidates in case a longer sequence could somehow score higher
            # (e.g., with positive gap/mismatch scores). Fall through to scoring.

        # Calculate alignment score using global alignment with custom scoring
        score = calculate_alignment_score(
            query_umi,
            candidate_str,
            mismatch_penalty,
            gap_penalty,
            match_score,
            min_score,
        )
        
        # Update best match if this scores higher and meets minimum threshold
        if score >= min_score and score > best_score:
            best_score = score
            best_match = candidate_str
            
            # If we found a perfect score for equal-length sequences, return immediately
            if score == max_possible_score and len(candidate_str) == len(query_umi):
                break
    
    return best_match


@lru_cache(maxsize=1000)
def calculate_alignment_score(
    seq1: str,
    seq2: str,
    mismatch_penalty: int,
    gap_penalty: int,
    match_score: int,
    min_score_cutoff: Optional[int] = None,
) -> int:
    """Calculate a global alignment score using simple dynamic programming.
    
    The scoring scheme is:
      - match: ``match_score``
      - mismatch: ``mismatch_penalty``
      - gap (insertion/deletion): ``gap_penalty``
    
    This function computes the Needleman–Wunsch global alignment score without
    traceback, returning only the optimal score.
    If ``min_score_cutoff`` is provided and the dynamic programming state proves
    the optimal score cannot reach the cutoff, the function returns a value
    strictly less than the cutoff early.
    """

    # Fast paths for empty inputs
    if not seq1 and not seq2:
        return 0
    if not seq1:
        score = len(seq2) * gap_penalty
        if min_score_cutoff is not None and score < min_score_cutoff:
            return min_score_cutoff - 1
        return score
    if not seq2:
        score = len(seq1) * gap_penalty
        if min_score_cutoff is not None and score < min_score_cutoff:
            return min_score_cutoff - 1
        return score

    len1 = len(seq1)
    len2 = len(seq2)

    # Initialize DP matrix with only two rows to reduce memory footprint
    prev_row = [0] * (len2 + 1)
    curr_row = [0] * (len2 + 1)

    # Base case: aligning empty prefix of seq1 with prefixes of seq2
    for j in range(1, len2 + 1):
        prev_row[j] = prev_row[j - 1] + gap_penalty

    # Early-abort feasibility checks apply only under standard scoring where
    # matches are non-negative and mismatches/gaps do not increase the score.
    can_prune = (
        min_score_cutoff is not None
        and match_score >= 0
        and mismatch_penalty <= 0
        and gap_penalty <= 0
    )

    # Global optimistic bound from the start
    if can_prune:
        max_initial = min(len1, len2) * match_score
        if max_initial < min_score_cutoff:  # impossible to reach cutoff
            return min_score_cutoff - 1

    for i in range(1, len1 + 1):
        # Base case for this row: aligning prefix of seq1 with empty prefix of seq2
        curr_row[0] = prev_row[0] + gap_penalty

        ch1 = seq1[i - 1]
        row_best = None
        for j in range(1, len2 + 1):
            ch2 = seq2[j - 1]
            diag = prev_row[j - 1] + (match_score if ch1 == ch2 else mismatch_penalty)
            up = prev_row[j] + gap_penalty      # gap in seq2 (deletion from seq1)
            left = curr_row[j - 1] + gap_penalty  # gap in seq1 (insertion into seq1)
            val = diag if diag >= up and diag >= left else (up if up >= left else left)
            curr_row[j] = val
            if row_best is None or val > row_best:
                row_best = val

        # Prune if even the most optimistic completion of this row can't hit cutoff
        if can_prune:
            # Remaining positions that could contribute positively as matches
            # For any position j in this completed row, the maximum extra is
            # (min(len1 - i, len2 - j)) * match_score. The maximum over j is
            # achieved at the smallest j (i.e., j=0), but j=0 may not be the
            # position of row_best. To keep this cheap and admissible, use the
            # most optimistic bound based on row_best at some j and allow up to
            # min(len1 - i, len2) future matches.
            optimistic_extra = min(len1 - i, len2) * match_score
            if (row_best or 0) + optimistic_extra < min_score_cutoff:
                return min_score_cutoff - 1

        # Prepare for next iteration
        prev_row, curr_row = curr_row, prev_row

    return prev_row[len2]


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
