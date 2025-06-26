"""Single UMI library statistics calculations.

This module provides classes for calculating various statistics on individual
UMI libraries including diversity metrics, efficiency measures, and error rates.
"""

import logging
from typing import Dict, List, Optional, Sequence, Any
import numpy as np
import pandas as pd
from ..umi import UMI
from .base import BaseStatistic
from .utils import split_counts_by_allowed_list

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


class UMIStats(BaseStatistic):
    """Base class for UMI statistics calculations.

    This class provides common functionality for UMI statistics calculations
    including count access and allowed list handling.
    """

    def __init__(
        self,
        umi: UMI,
        use_corrected: bool = True,
        allowed_list: Optional[List[str]] = None,
        **kwargs: Any,
    ) -> None:
        """Initialize UMI stats calculator.

        Parameters
        ----------
        umi : UMI
            UMI object to calculate statistics on
        use_corrected : bool, default=True
            If True, use corrected counts. If False, use original counts.
        allowed_list : Optional[List[str]], default=None
            Optional list of allowed UMIs for filtering
        **kwargs : Any
            Additional arguments passed to parent class
        """
        super().__init__(umi, **kwargs)
        self.use_corrected = use_corrected
        self.allowed_list = allowed_list

    @property
    def _allowed_list(self) -> Optional[List[bytes]]:
        """Get allowed list in bytes format.

        Returns
        -------
        Optional[List[bytes]]
            Allowed list converted to bytes, or None if not provided
        """
        if self.allowed_list:
            return [umi.encode("ascii") for umi in self.allowed_list]
        return None

    @property
    def _counts(self) -> Dict[bytes, int]:
        """Get counts based on use_corrected setting.

        Returns
        -------
        Dict[bytes, int]
            Either corrected or original counts
        """
        if self.use_corrected:
            return self.umi.corrected_counts
        return self.umi._counts


class GiniCoefficient(UMIStats):
    """Calculate Gini coefficient for UMI counts.

    The Gini coefficient measures inequality in the distribution of UMI counts.
    A value of 0 indicates perfect equality, while a value of 1 indicates
    maximum inequality.
    """

    def __init__(
        self,
        umi: UMI,
        use_corrected: bool = True,
        allowed_list: Optional[List[str]] = None,
    ) -> None:
        """Initialize Gini coefficient calculator.

        Parameters
        ----------
        umi : UMI
            UMI object to calculate statistics on
        use_corrected : bool, default=True
            If True, use corrected counts. If False, use original counts.
        allowed_list : Optional[List[str]], default=None
            Optional list of allowed keys. If provided, missing keys will be
            treated as having zero counts in the calculation.
        """
        super().__init__(umi, use_corrected, allowed_list)

    @staticmethod
    def calculate_gini(counts: Sequence[float]) -> Optional[float]:
        """Calculate Gini coefficient from a sequence of counts.

        Parameters
        ----------
        counts : Sequence[float]
            Sequence of count values

        Returns
        -------
        Optional[float]
            Gini coefficient or None if calculation is not possible

        Notes
        -----
        The Gini coefficient is calculated using the formula:
        G = (2 * sum(i * y_i)) / (n * sum(y_i)) - (n + 1) / n
        where i is the rank and y_i is the count value
        """
        if not counts:
            return None

        # Sort counts in ascending order
        sorted_counts = sorted(counts)
        n = len(sorted_counts)

        # If all counts are zero, return None
        total = sum(sorted_counts)
        if total == 0:
            return None

        # Calculate the Lorenz curve
        index = list(range(1, n + 1))
        gini = (
            (2 * sum(i * y for i, y in zip(index, sorted_counts))) / (n * total)
        ) - ((n + 1) / n)

        return gini

    def run(self) -> Optional[float]:
        """Calculate the Gini coefficient for the UMI object.

        Returns
        -------
        Optional[float]
            Gini coefficient or None if calculation is not possible
        """
        if self.allowed_list:
            counts, _, _ = split_counts_by_allowed_list(
                self._counts, self._allowed_list, add_missing=True
            )
        else:
            counts = self._counts

        result = GiniCoefficient.calculate_gini(list(counts.values()))
        logger.debug(f"Calculated Gini coefficient: {result}")
        return result


class ShannonDiversity(UMIStats):
    """Calculate Shannon diversity index for UMI counts.

    The Shannon diversity index measures the diversity and evenness of UMI
    distribution. Higher values indicate greater diversity.
    """

    def __init__(
        self,
        umi: UMI,
        use_corrected: bool = True,
        base: float = 2.0,
        allowed_list: Optional[List[str]] = None,
    ) -> None:
        """Initialize Shannon diversity calculator.

        Parameters
        ----------
        umi : UMI
            UMI object to calculate statistics on
        use_corrected : bool, default=True
            If True, use corrected counts. If False, use original counts.
        base : float, default=2.0
            Base of the logarithm (default: 2.0 for bits)
        allowed_list : Optional[List[str]], default=None
            Optional list of allowed UMIs for filtering
        """
        super().__init__(umi, use_corrected, allowed_list)
        self.base = base

    @staticmethod
    def calculate_shannon(
        counts: Sequence[float], base: float = 2.0
    ) -> Optional[float]:
        """Calculate Shannon diversity index from a sequence of counts.

        Parameters
        ----------
        counts : Sequence[float]
            Sequence of count values
        base : float, default=2.0
            Base of the logarithm (default: 2.0 for bits)

        Returns
        -------
        Optional[float]
            Shannon diversity index or None if calculation is not possible

        Notes
        -----
        The Shannon diversity index is calculated using the formula:
        H = -sum(p_i * log(p_i)) / log(base)
        where p_i is the proportion of each count
        """
        if not counts:
            return None

        # Convert to numpy array for efficient calculation
        counts_array = np.array([c for c in counts if c > 0])
        total = np.sum(counts_array)

        if total == 0:
            return None

        # Calculate proportions
        proportions = counts_array / total

        # Calculate Shannon diversity
        shannon = -np.sum(proportions * np.log(proportions) / np.log(base))

        return shannon

    def run(self) -> Optional[float]:
        """Calculate the Shannon diversity index for the UMI object.

        Returns
        -------
        Optional[float]
            Shannon diversity index or None if calculation is not possible
        """
        if self.allowed_list:
            counts, _, _ = split_counts_by_allowed_list(
                self._counts, self._allowed_list, add_missing=False
            )
        else:
            counts = self._counts

        result = ShannonDiversity.calculate_shannon(list(counts.values()), self.base)
        logger.debug(f"Calculated Shannon diversity (base {self.base}): {result}")
        return result


class SimpsonDiversity(UMIStats):
    """Calculate Simpson's diversity index for UMI counts.

    Simpson's diversity index measures the probability that two randomly
    selected UMIs are different. Higher values indicate greater diversity.
    """

    def __init__(
        self,
        umi: UMI,
        use_corrected: bool = True,
        allowed_list: Optional[List[str]] = None,
    ) -> None:
        """Initialize Simpson's diversity calculator.

        Parameters
        ----------
        umi : UMI
            UMI object to calculate statistics on
        use_corrected : bool, default=True
            If True, use corrected counts. If False, use original counts.
        allowed_list : Optional[List[str]], default=None
            Optional list of allowed keys. If provided, missing keys will be
            treated as having zero counts in the calculation.
        """
        super().__init__(umi, use_corrected, allowed_list)

    @staticmethod
    def calculate_simpson(counts: Sequence[float]) -> Optional[float]:
        """Calculate Simpson's diversity index from a sequence of counts.

        Parameters
        ----------
        counts : Sequence[float]
            Sequence of count values

        Returns
        -------
        Optional[float]
            Simpson's diversity index or None if calculation is not possible

        Notes
        -----
        Simpson's diversity index is calculated using the formula:
        D = 1 - sum(p_i^2)
        where p_i is the proportion of each count
        """
        if not counts:
            return None

        # Convert dict_values to list before numpy array conversion
        counts_array = np.array(list(counts))
        total = np.sum(counts_array)

        if total == 0:
            return None

        # Calculate proportions
        proportions = counts_array / total

        # Calculate Simpson's diversity (1 - D)
        simpson = 1 - np.sum(proportions**2)

        return simpson

    def run(self) -> Optional[float]:
        """Calculate Simpson's diversity index for the UMI object.

        Returns
        -------
        Optional[float]
            Simpson's diversity index or None if calculation is not possible
        """
        if self.allowed_list:
            counts, _, _ = split_counts_by_allowed_list(
                self._counts, self._allowed_list, add_missing=False
            )
        else:
            counts = self._counts

        result = SimpsonDiversity.calculate_simpson(list(counts.values()))
        logger.debug(f"Calculated Simpson diversity: {result}")
        return result


class UMIRecoveryRate(UMIStats):
    """Calculate UMI recovery rate.

    The UMI recovery rate measures the fraction of expected UMIs that were
    actually observed in the data.
    """

    def __init__(
        self,
        umi: UMI,
        use_corrected: bool = True,
        allowed_list: Optional[List[str]] = None,
    ) -> None:
        """Initialize UMI recovery rate calculator.

        Parameters
        ----------
        umi : UMI
            UMI object to calculate statistics on
        use_corrected : bool, default=True
            If True, use corrected counts. If False, use original counts.
        allowed_list : Optional[List[str]], default=None
            Optional list of allowed UMIs. If provided, recovery rate is
            calculated as the ratio of observed allowed UMIs to total allowed UMIs.
            If not provided, assumes exponential distribution and calculates
            theoretical recovery rate.
        """
        super().__init__(umi, use_corrected, allowed_list)

    @staticmethod
    def calculate_recovery_rate_limited(
        counts: Dict[bytes, int], allowed_list: List[bytes]
    ) -> Optional[float]:
        """Calculate UMI recovery rate from observed and total unique counts.

        Parameters
        ----------
        counts : Dict[bytes, int]
            Dictionary of UMIs and their counts
        allowed_list : List[bytes]
            List of allowed UMIs

        Returns
        -------
        Optional[float]
            Recovery rate as a fraction between 0 and 1
        """
        _, _, missing = split_counts_by_allowed_list(
            counts, allowed_list, add_missing=True
        )
        return (len(allowed_list) - len(missing)) / len(allowed_list)

    def run(self) -> Optional[float]:
        """Calculate UMI recovery rate for the UMI object.

        Returns
        -------
        Optional[float]
            UMI recovery rate or None if calculation is not possible
        """
        if not self.allowed_list:
            return None

        result = UMIRecoveryRate.calculate_recovery_rate_limited(
            self._counts, self._allowed_list
        )
        logger.debug(f"Calculated UMI recovery rate: {result}")
        return result


class UMIEfficiencyRate(UMIStats):
    """Calculate UMI efficiency rate (fraction of reads contributing to allowed UMIs).

    The UMI efficiency rate measures the fraction of total reads that contributed
    to allowed UMIs rather than being wasted on unwanted UMIs.
    """

    def __init__(
        self,
        umi: UMI,
        use_corrected: bool = True,
        allowed_list: Optional[List[str]] = None,
    ) -> None:
        """Initialize UMI efficiency rate calculator.

        Parameters
        ----------
        umi : UMI
            UMI object to calculate statistics on
        allowed_list : List[str]
            List of allowed UMIs. Efficiency rate is calculated as the
            fraction of reads that contributed to these allowed UMIs.
        use_corrected : bool, default=True
            If True, use corrected counts. If False, use original counts.
        """
        super().__init__(umi, use_corrected, allowed_list)

    @staticmethod
    def calculate_efficiency_rate(
        counts: Dict[bytes, int], allowed_list: List[bytes]
    ) -> Optional[float]:
        """Calculate UMI efficiency rate from allowed and total reads.

        Parameters
        ----------
        counts : Dict[bytes, int]
            Dictionary of UMIs and their counts
        allowed_list : List[bytes]
            List of allowed UMIs

        Returns
        -------
        Optional[float]
            UMI efficiency rate or None if calculation is not possible

        Notes
        -----
        UMI efficiency rate is calculated as allowed_reads / total_reads
        """
        if not counts:
            return None

        counts, banned, _ = split_counts_by_allowed_list(
            counts, allowed_list, add_missing=True
        )

        wanted_reads = sum(counts.values())
        total_reads = sum(counts.values()) + sum(banned.values())

        if total_reads == 0:
            return None

        return wanted_reads / total_reads

    def run(self) -> Optional[float]:
        """Calculate UMI efficiency rate for the UMI object.

        Returns
        -------
        Optional[float]
            UMI efficiency rate or None if calculation is not possible
        """
        if not self.allowed_list:
            return None

        result = UMIEfficiencyRate.calculate_efficiency_rate(
            self._counts, self._allowed_list
        )
        logger.debug(f"Calculated UMI efficiency rate: {result}")
        return result


class UMIErrorRate(BaseStatistic):
    """Calculate UMI error rate based on mismatches between original and corrected UMIs.

    The UMI error rate measures the average number of mismatches per read
    between original and corrected UMI sequences.
    """

    def __init__(self, umi: UMI, **kwargs: Any) -> None:
        """Initialize UMI error rate calculator.

        Parameters
        ----------
        umi : UMI
            UMI object to calculate statistics on
        **kwargs : Any
            Additional arguments passed to parent class
        """
        super().__init__(umi, **kwargs)

    @staticmethod
    def hamming_distance(seq1: bytes, seq2: bytes) -> int:
        """Calculate Hamming distance between two sequences.

        Parameters
        ----------
        seq1 : bytes
            First sequence
        seq2 : bytes
            Second sequence

        Returns
        -------
        int
            Number of mismatches between sequences

        Raises
        ------
        ValueError
            If sequences are of different lengths
        """
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must be of the same length")
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

    @staticmethod
    def calculate_error_rate(
        mapping: Dict[bytes, bytes], counts: Dict[bytes, int]
    ) -> Optional[float]:
        """Calculate UMI error rate from mismatches and total reads.

        Parameters
        ----------
        mapping : Dict[bytes, bytes]
            Dictionary of original and corrected UMIs
        counts : Dict[bytes, int]
            Dictionary of UMIs and their counts

        Returns
        -------
        Optional[float]
            UMI error rate or None if calculation is not possible

        Notes
        -----
        UMI error rate is calculated as total mismatches / total reads
        """
        # Calculate total mismatches
        if not mapping or not counts:
            return None

        total_mismatches = 0
        for original, corrected in mapping.items():
            if original != corrected:  # Only calculate mismatches if UMI was corrected
                mismatches = UMIErrorRate.hamming_distance(original, corrected)
                total_mismatches += mismatches * counts[original]

        total_reads = sum(counts.values())

        return total_mismatches / total_reads

    def run(self) -> Optional[float]:
        """Calculate UMI error rate for the UMI object.

        Returns
        -------
        Optional[float]
            UMI error rate or None if calculation is not possible
        """
        if not self.umi._mapping or not self.umi._counts:
            return None

        result = UMIErrorRate.calculate_error_rate(self.umi._mapping, self.umi._counts)
        logger.debug(f"Calculated UMI error rate: {result}")
        return result


class UMIRedundancy(UMIStats):
    """Calculate UMI redundancy (average reads per unique UMI).

    UMI redundancy measures the average number of reads per unique UMI,
    indicating how much sequencing depth is available per unique sequence.
    """

    def __init__(
        self,
        umi: UMI,
        use_corrected: bool = True,
        allowed_list: Optional[List[str]] = None,
    ) -> None:
        """Initialize UMI redundancy calculator.

        Parameters
        ----------
        umi : UMI
            UMI object to calculate statistics on
        use_corrected : bool, default=True
            If True, use corrected counts. If False, use original counts.
        allowed_list : Optional[List[str]], default=None
            Optional list of allowed UMIs for filtering
        """
        super().__init__(umi, use_corrected, allowed_list)

    @staticmethod
    def calculate_redundancy(counts: Dict[bytes, int]) -> Optional[float]:
        """Calculate UMI redundancy from total and unique counts.

        Parameters
        ----------
        counts : Dict[bytes, int]
            Dictionary of UMIs and their counts

        Returns
        -------
        Optional[float]
            UMI redundancy or None if calculation is not possible

        Notes
        -----
        UMI redundancy is calculated as total_reads / unique_umis
        """
        total_reads = sum(counts.values())
        unique_umis = len(counts)

        if unique_umis == 0:
            return None

        return total_reads / unique_umis

    def run(self) -> Optional[float]:
        """Calculate UMI redundancy for the UMI object.

        Returns
        -------
        Optional[float]
            UMI redundancy or None if calculation is not possible
        """
        if self.allowed_list:
            counts, _, _ = split_counts_by_allowed_list(
                self._counts, self._allowed_list, add_missing=False
            )
        else:
            counts = self._counts

        if not counts:
            return None

        result = UMIRedundancy.calculate_redundancy(counts)
        logger.debug(f"Calculated UMI redundancy: {result}")
        return result


# Copyright (C) 2025, SC Barrera, Drs DVK & WND. All Rights Reserved.
