"""Nearest neighbor matching functionality for UMIs.

This module provides functions for finding the closest matching UMI from a list
of allowed UMIs using a configurable global alignment scoring system.
"""

import logging
from typing import Optional, Tuple, List, Iterable, Dict, Set
from functools import lru_cache

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


class NearestUMIFinder:
    """Class-based nearest neighbor finder for UMIs.

    Holds the allowed UMI list and alignment scoring parameters.
    Provides methods to find the closest UMI(s) for a query or a batch of queries.
    """

    def __init__(
        self,
        allowed_list: Optional[Iterable[str]] = None,
        mismatch_penalty: int = -1,
        gap_penalty: int = -3,
        match_score: int = 1,
        min_score: int = 0,
        use_prescreen: bool = True,
        kmer_size: int = 3,
        min_kmer_overlap: int = 1,
    ) -> None:
        self._allowed_umis_tuple: Tuple[str, ...] = tuple(allowed_list or ())
        self._allowed_set: Set[str] = set(self._allowed_umis_tuple)
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty
        self.match_score = match_score
        self.min_score = min_score
        self._use_prescreen = use_prescreen
        self._kmer_size = kmer_size
        self._min_kmer_overlap = min_kmer_overlap

        # Pre-computed k-mer sets and inverted index for prescreening
        self._kmer_sets: Dict[str, frozenset] = {}
        if self._use_prescreen:
            self._rebuild_kmer_index()

    @property
    def allowed_list(self) -> List[str]:
        return list(self._allowed_umis_tuple)

    def set_allowed_list(self, allowed_list: Iterable[str]) -> None:
        """Replace the allowed UMI list."""
        self._allowed_umis_tuple = tuple(allowed_list or ())
        self._allowed_set = set(self._allowed_umis_tuple)
        if self._use_prescreen:
            self._rebuild_kmer_index()

    def add_allowed(self, umi: str) -> None:
        """Add a single UMI to the allowed list."""
        if umi in self._allowed_umis_tuple:
            return
        self._allowed_umis_tuple = tuple(list(self._allowed_umis_tuple) + [umi])
        self._allowed_set.add(umi)
        if self._use_prescreen:
            self._kmer_sets[umi] = self._kmerize(umi)

    def find(self, query_umi: str) -> Optional[List[str]]:
        """Find closest allowed UMI(s) for a single query.

        Returns a list of best-scoring matches or None if below ``min_score``.
        """
        # Validate inputs
        if not self._allowed_umis_tuple or not query_umi:
            return None

        # Exact match short-circuit using set membership
        if query_umi in self._allowed_set:
            return [query_umi]

        # Initialize tracking of best matches
        best_matches: List[str] = []
        best_score: int = self.min_score - 1

        # Determine candidate subset via k-mer prescreen if enabled
        if self._use_prescreen and self._kmer_size and self._kmer_size > 0:
            kmers = self._kmerize(query_umi)
        else:
            kmers = None

        # Iterate candidates
        for candidate_str in self._iter_candidates(query_kmers=kmers):
            # Exact match short-circuit
            if self._is_exact_match(candidate_str, query_umi):
                return [candidate_str]

            score = self._score_candidate(query_umi, candidate_str, best_score)
            if self._below_threshold(score):
                continue

            if score == best_score:
                best_matches.append(candidate_str)
            elif score > best_score:
                best_score = score
                best_matches = [candidate_str]

        return best_matches or None

    def find_many(self, query_umis: Iterable[str]) -> List[Optional[List[str]]]:
        """Find closest allowed UMI(s) for multiple queries in order."""
        return [self.find(q) for q in query_umis]

    # ---------- Private helpers ----------
    def _iter_candidates(self, query_kmers: Optional[frozenset] = None) -> Iterable[str]:
        if query_kmers:
            for cand, cand_kmers in self._kmer_sets.items():
                if len(cand_kmers.intersection(query_kmers)) >= self._min_kmer_overlap:
                    yield cand
        else:
            yield from self._allowed_umis_tuple

    def _is_exact_match(self, candidate: str, query: str) -> bool:
        return candidate == query

    def _score_candidate(self, query: str, candidate: str, current_best: int) -> int:
        return NearestUMIFinder._calculate_alignment_score(
            query,
            candidate,
            self.mismatch_penalty,
            self.gap_penalty,
            self.match_score,
            current_best,
        )

    def _below_threshold(self, score: int) -> bool:
        return score < self.min_score

    def _kmerize(self, seq: str) -> frozenset:
        """Return a frozenset of sliding k-mers for a sequence.

        If the sequence is shorter than k, returns a singleton set of the sequence
        to preserve some overlap semantics.
        """
        k = self._kmer_size
        if k <= 0:
            return frozenset()
        if len(seq) < k:
            return frozenset([seq])
        return frozenset(seq[i : i + k] for i in range(len(seq) - k + 1))

    def _rebuild_kmer_index(self) -> None:
        """Rebuild k-mer sets and inverted index from the allowed list."""
        self._kmer_sets = {}
        if not self._allowed_umis_tuple:
            return
        for cand in self._allowed_umis_tuple:
            self._kmer_sets[cand] = self._kmerize(cand)

    @staticmethod
    def _calculate_alignment_score(
        seq1: str,
        seq2: str,
        mismatch_penalty: int,
        gap_penalty: int,
        match_score: int,
        min_score_cutoff: Optional[int] = None,
    ) -> int:
        """Alignment score implementation (Needleman–Wunsch without traceback).

        Mirrors the behavior of the module-level function.
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

        # Early-abort feasibility checks apply only under standard scoring
        can_prune = (
            min_score_cutoff is not None
            and match_score >= 0
            and mismatch_penalty <= 0
            and gap_penalty <= 0
        )

        # Global optimistic bound from the start
        if can_prune:
            max_initial = min(len1, len2) * match_score
            if max_initial < min_score_cutoff:
                return min_score_cutoff - 1

        for i in range(1, len1 + 1):
            # Base case for this row: aligning prefix of seq1 with empty prefix of seq2
            curr_row[0] = prev_row[0] + gap_penalty

            ch1 = seq1[i - 1]
            row_best = None
            for j in range(1, len2 + 1):
                ch2 = seq2[j - 1]
                diag = prev_row[j - 1] + (match_score if ch1 == ch2 else mismatch_penalty)
                up = prev_row[j] + gap_penalty
                left = curr_row[j - 1] + gap_penalty
                val = diag if diag >= up and diag >= left else (up if up >= left else left)
                curr_row[j] = val
                if row_best is None or val > row_best:
                    row_best = val

            # Prune if even the most optimistic completion of this row can't hit cutoff
            if can_prune:
                optimistic_extra = min(len1 - i, len2) * match_score
                if (row_best or 0) + optimistic_extra < min_score_cutoff:
                    return min_score_cutoff - 1

            # Prepare for next iteration
            prev_row, curr_row = curr_row, prev_row

        return prev_row[len2]


@lru_cache(maxsize=50_000)
def find_closest_umi(
    allowed_umis: Tuple[str, ...], 
    query_umi: str, 
    mismatch_penalty: int = -1,
    gap_penalty: int = -3,
    match_score: int = 1,
    min_score: int = 0
) -> List[str]:
    """Wrapper that delegates to ``NearestUMIFinder`` for nearest UMI search."""
    finder = NearestUMIFinder(
        allowed_list=allowed_umis,
        mismatch_penalty=mismatch_penalty,
        gap_penalty=gap_penalty,
        match_score=match_score,
        min_score=min_score,
    )
    return finder.find(query_umi)


@lru_cache(maxsize=50_000)
def calculate_alignment_score(
    seq1: str,
    seq2: str,
    mismatch_penalty: int,
    gap_penalty: int,
    match_score: int,
    min_score_cutoff: Optional[int] = None,
) -> int:
    """Wrapper delegating to class implementation (cached)."""
    return NearestUMIFinder._calculate_alignment_score(
        seq1,
        seq2,
        mismatch_penalty,
        gap_penalty,
        match_score,
        min_score_cutoff,
    )


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
