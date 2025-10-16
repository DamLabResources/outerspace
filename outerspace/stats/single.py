"""Single UMI library statistics calculations.

This module provides classes for calculating various statistics on individual
UMI libraries including diversity metrics, efficiency measures, and error rates.
"""

import csv
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Any, Set
import numpy as np
import pandas as pd
from ..umi import UMI
from .base import BaseStatistic
from .utils import split_counts_by_allowed_list

# Set up logging
logger = logging.getLogger(__name__)

# Increase CSV field size limit to handle large fields
csv.field_size_limit(sys.maxsize)

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


def _read_allowed_list(filepath: str) -> List[str]:
    """Read allowed values from a text file.

    Parameters
    ----------
    filepath : str
        Path to text file containing allowed values

    Returns
    -------
    List[str]
        List of allowed values with empty lines filtered out
    """
    with open(filepath, "r") as f:
        return [line.strip() for line in f if line.strip()]


def _aggregate_counts_from_csv(
    input_file: str,
    key_column: str,
    barcode_column: Optional[str] = None,
    sep: str = ","
) -> Dict[str, int]:
    """Aggregate counts from collapse output CSV.

    This function reads a collapse output CSV and aggregates the data to produce
    counts suitable for creating a UMI object for statistics calculation.

    Parameters
    ----------
    input_file : str
        Path to input CSV file (collapse output)
    key_column : str
        Column containing the keys (e.g., corrected sequence)
    barcode_column : Optional[str], default=None
        If provided, count unique values in this column per key (e.g., unique UMIs per sequence).
        If None, count occurrences of each key (e.g., number of reads per sequence).
    sep : str, default=','
        CSV separator

    Returns
    -------
    Dict[str, int]
        Dictionary mapping keys to their counts

    Raises
    ------
    ValueError
        If required columns are not found in the CSV file
    """
    if barcode_column:
        # Count unique barcodes per key
        key_barcodes: Dict[str, Set[str]] = defaultdict(set)
        
        with open(input_file, "r") as f:
            reader = csv.DictReader(f, delimiter=sep)
            
            if key_column not in reader.fieldnames:
                raise ValueError(f"Column '{key_column}' not found in {input_file}")
            if barcode_column not in reader.fieldnames:
                raise ValueError(f"Column '{barcode_column}' not found in {input_file}")
            
            for row in reader:
                key = row[key_column]
                barcode = row[barcode_column]
                
                if key and barcode:  # Skip empty values
                    key_barcodes[key].add(barcode)
        
        # Convert to counts
        counts = {key: len(barcodes) for key, barcodes in key_barcodes.items()}
        logger.debug(
            f"Aggregated {sum(counts.values())} unique barcodes across {len(counts)} keys"
        )
        return counts
    else:
        # Count occurrences of each key
        key_counts: Dict[str, int] = defaultdict(int)
        
        with open(input_file, "r") as f:
            reader = csv.DictReader(f, delimiter=sep)
            
            if key_column not in reader.fieldnames:
                raise ValueError(f"Column '{key_column}' not found in {input_file}")
            
            for row in reader:
                key = row[key_column]
                
                if key:  # Skip empty values
                    key_counts[key] += 1
        
        counts = dict(key_counts)
        logger.debug(
            f"Counted {sum(counts.values())} total occurrences across {len(counts)} keys"
        )
        return counts


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


    # TODO: Remove allowed_list logic across all stats classes, it is handled elsewhere
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

    @classmethod
    def _from_step(cls, input_file: str, sep: str = ",", **step_params: Any) -> Optional[float]:
        """Create statistic from step parameters.

        Parameters
        ----------
        input_file : str
            Path to collapse output CSV file
        sep : str, default=','
            CSV separator
        **step_params : Any
            Should include:
            - key_column: str - Column containing keys (required)
            - barcode_column: str (optional) - If provided, count unique barcodes per key
            - allowed_list: str (optional) - Path to allowed list file
            - use_corrected: bool (optional) - Whether to use corrected counts

        Returns
        -------
        Optional[float]
            Calculated Gini coefficient
        """
        key_column = step_params.get("key_column")
        barcode_column = step_params.get("barcode_column")
        allowed_list_path = step_params.get("allowed_list")
        use_corrected = step_params.get("use_corrected", False)  # No correction needed, data is pre-collapsed

        if not key_column:
            raise ValueError("key_column is required for GiniCoefficient")

        # Load allowed list if provided
        allowed_list = None
        if allowed_list_path:
            allowed_list = _read_allowed_list(allowed_list_path)

        # Aggregate counts from collapse output
        counts_dict = _aggregate_counts_from_csv(input_file, key_column, barcode_column, sep)
        
        # Create UMI object from aggregated counts
        keys = list(counts_dict.keys())
        counts = list(counts_dict.values())
        umi = UMI()
        for key, count in zip(keys, counts):
            umi.consume(key, count)
        # Set mapping for pre-collapsed data (each key maps to itself)
        umi._mapping = {k.encode() if isinstance(k, str) else k: k.encode() if isinstance(k, str) else k for k in keys}
        umi._corrected_counts = umi._counts.copy()

        return cls.calculate(umi, use_corrected=use_corrected, allowed_list=allowed_list)


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

    @classmethod
    def _from_step(cls, input_file: str, sep: str = ",", **step_params: Any) -> Optional[float]:
        """Create statistic from step parameters.

        Parameters
        ----------
        input_file : str
            Path to collapse output CSV file
        sep : str, default=','
            CSV separator
        **step_params : Any
            Should include:
            - key_column: str - Column containing keys (required)
            - barcode_column: str (optional) - If provided, count unique barcodes per key
            - allowed_list: str (optional) - Path to allowed list file
            - use_corrected: bool (optional) - Whether to use corrected counts
            - base: float (optional) - Base for logarithm (default 2.0)

        Returns
        -------
        Optional[float]
            Calculated Shannon diversity
        """
        key_column = step_params.get("key_column")
        barcode_column = step_params.get("barcode_column")
        allowed_list_path = step_params.get("allowed_list")
        use_corrected = step_params.get("use_corrected", False)
        base = step_params.get("base", 2.0)

        if not key_column:
            raise ValueError("key_column is required for ShannonDiversity")

        # Load allowed list if provided
        allowed_list = None
        if allowed_list_path:
            allowed_list = _read_allowed_list(allowed_list_path)

        # Aggregate counts from collapse output
        counts_dict = _aggregate_counts_from_csv(input_file, key_column, barcode_column, sep)
        
        # Create UMI object from aggregated counts
        keys = list(counts_dict.keys())
        counts = list(counts_dict.values())
        umi = UMI()
        for key, count in zip(keys, counts):
            umi.consume(key, count)
        umi._mapping = {k.encode() if isinstance(k, str) else k: k.encode() if isinstance(k, str) else k for k in keys}
        umi._corrected_counts = umi._counts.copy()

        return cls.calculate(umi, use_corrected=use_corrected, allowed_list=allowed_list, base=base)


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

    @classmethod
    def _from_step(cls, input_file: str, sep: str = ",", **step_params: Any) -> Optional[float]:
        """Create statistic from step parameters.

        Parameters
        ----------
        input_file : str
            Path to collapse output CSV file
        sep : str, default=','
            CSV separator
        **step_params : Any
            Should include:
            - key_column: str - Column containing keys (required)
            - barcode_column: str (optional) - If provided, count unique barcodes per key
            - allowed_list: str (optional) - Path to allowed list file
            - use_corrected: bool (optional) - Whether to use corrected counts

        Returns
        -------
        Optional[float]
            Calculated Simpson diversity
        """
        key_column = step_params.get("key_column")
        barcode_column = step_params.get("barcode_column")
        allowed_list_path = step_params.get("allowed_list")
        use_corrected = step_params.get("use_corrected", False)

        if not key_column:
            raise ValueError("key_column is required for SimpsonDiversity")

        # Load allowed list if provided
        allowed_list = None
        if allowed_list_path:
            allowed_list = _read_allowed_list(allowed_list_path)

        # Aggregate counts from collapse output
        counts_dict = _aggregate_counts_from_csv(input_file, key_column, barcode_column, sep)
        
        # Create UMI object from aggregated counts
        keys = list(counts_dict.keys())
        counts = list(counts_dict.values())
        umi = UMI()
        for key, count in zip(keys, counts):
            umi.consume(key, count)
        umi._mapping = {k.encode() if isinstance(k, str) else k: k.encode() if isinstance(k, str) else k for k in keys}
        umi._corrected_counts = umi._counts.copy()

        return cls.calculate(umi, use_corrected=use_corrected, allowed_list=allowed_list)

# TODO: Add a class to calculcate the Hill number



class UMIRecoveryRate(UMIStats):
    """Calculate UMI recovery rate.

    The UMI recovery rate measures the fraction of expected UMIs that were
    actually observed in the data.
    """

    def __init__(
        self,
        umi: UMI,
        use_corrected: bool = True,
        allowed_list: Optional[List[str]] = None, # This allowed list should stay
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

    @classmethod
    def _from_step(cls, input_file: str, sep: str = ",", **step_params: Any) -> Optional[float]:
        """Create statistic from step parameters.

        Parameters
        ----------
        input_file : str
            Path to input CSV file
        sep : str, default=','
            CSV separator
        **step_params : Any
            Should include:
            - key_column: str - Column containing keys (required)
            - barcode_column: str (optional) - If provided, count unique barcodes per key
            - allowed_list: str - Path to allowed list file (required)
            - use_corrected: bool (optional) - Whether to use corrected counts

        Returns
        -------
        Optional[float]
            Calculated UMI recovery rate
        """
        key_column = step_params.get("key_column")
        barcode_column = step_params.get("barcode_column")
        allowed_list_path = step_params.get("allowed_list")
        use_corrected = step_params.get("use_corrected", False)

        if not key_column:
            raise ValueError("key_column is required for UMIRecoveryRate")
        if not allowed_list_path:
            raise ValueError("allowed_list is required for UMIRecoveryRate")

        # Load allowed list
        allowed_list = _read_allowed_list(allowed_list_path)

        # Aggregate counts from collapse output
        counts_dict = _aggregate_counts_from_csv(input_file, key_column, barcode_column, sep)
        
        # Create UMI object from aggregated counts
        keys = list(counts_dict.keys())
        counts = list(counts_dict.values())
        umi = UMI()
        for key, count in zip(keys, counts):
            umi.consume(key, count)
        umi._mapping = {k.encode() if isinstance(k, str) else k: k.encode() if isinstance(k, str) else k for k in keys}
        umi._corrected_counts = umi._counts.copy()

        return cls.calculate(umi, use_corrected=use_corrected, allowed_list=allowed_list)


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

    @classmethod
    def _from_step(cls, input_file: str, sep: str = ",", **step_params: Any) -> Optional[float]:
        """Create statistic from step parameters.

        Parameters
        ----------
        input_file : str
            Path to input CSV file
        sep : str, default=','
            CSV separator
        **step_params : Any
            Should include:
            - key_column: str - Column containing keys (required)
            - barcode_column: str (optional) - If provided, count unique barcodes per key
            - allowed_list: str - Path to allowed list file (required)
            - use_corrected: bool (optional) - Whether to use corrected counts

        Returns
        -------
        Optional[float]
            Calculated UMI efficiency rate
        """
        key_column = step_params.get("key_column")
        barcode_column = step_params.get("barcode_column")
        allowed_list_path = step_params.get("allowed_list")
        use_corrected = step_params.get("use_corrected", False)

        if not key_column:
            raise ValueError("key_column is required for UMIEfficiencyRate")
        if not allowed_list_path:
            raise ValueError("allowed_list is required for UMIEfficiencyRate")

        # Load allowed list
        allowed_list = _read_allowed_list(allowed_list_path)

        # Aggregate counts from collapse output
        counts_dict = _aggregate_counts_from_csv(input_file, key_column, barcode_column, sep)
        
        # Create UMI object from aggregated counts
        keys = list(counts_dict.keys())
        counts = list(counts_dict.values())
        umi = UMI()
        for key, count in zip(keys, counts):
            umi.consume(key, count)
        umi._mapping = {k.encode() if isinstance(k, str) else k: k.encode() if isinstance(k, str) else k for k in keys}
        umi._corrected_counts = umi._counts.copy()

        return cls.calculate(umi, use_corrected=use_corrected, allowed_list=allowed_list)

# TODO: Redundant, remove this class
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

# TODO: Adjust this to calculate the error rate based on the score from nearest.NearestUMIFinder._calculate_alignment_score
class ErrorRate(BaseStatistic):
    """Calculate error rate by comparing original and corrected columns.

    This class compares two columns from a CSV file to calculate the
    error rate based on Hamming distance between corresponding values.
    """

    def __init__(self, umi: UMI, **kwargs: Any) -> None:
        """Initialize error rate calculator.

        Parameters
        ----------
        umi : UMI
            UMI object (not used in this implementation)
        **kwargs : Any
            Additional arguments (not used)
        """
        super().__init__(umi, **kwargs)

    def run(self) -> Optional[float]:
        """Not implemented for this class.

        Use _from_step classmethod instead.
        """
        raise NotImplementedError("Use ErrorRate._from_step() instead")

    @staticmethod
    def hamming_distance(seq1: str, seq2: str) -> int:
        """Calculate Hamming distance between two sequences.

        Parameters
        ----------
        seq1 : str
            First sequence
        seq2 : str
            Second sequence

        Returns
        -------
        int
            Number of mismatches between sequences
        """
        if len(seq1) != len(seq2):
            raise ValueError("Sequences must be of the same length")
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

    @classmethod
    def _from_step(cls, input_file: str, sep: str = ",", **step_params: Any) -> Optional[float]:
        """Create statistic from step parameters.

        Parameters
        ----------
        input_file : str
            Path to input CSV file
        sep : str, default=','
            CSV separator
        **step_params : Any
            Should include:
            - original_column: str - Column with original values
            - corrected_column: str - Column with corrected values

        Returns
        -------
        Optional[float]
            Calculated error rate (errors per position)
        """
        original_column = step_params.get("original_column")
        corrected_column = step_params.get("corrected_column")

        if not original_column or not corrected_column:
            raise ValueError("Both original_column and corrected_column are required for ErrorRate")

        total_mismatches = 0
        total_positions = 0

        with open(input_file, "r") as f:
            reader = csv.DictReader(f, delimiter=sep)
            
            if original_column not in reader.fieldnames:
                raise ValueError(f"Column '{original_column}' not found in {input_file}")
            if corrected_column not in reader.fieldnames:
                raise ValueError(f"Column '{corrected_column}' not found in {input_file}")

            for row in reader:
                original = row[original_column]
                corrected = row[corrected_column]
                
                if not original or not corrected:
                    continue

                if len(original) != len(corrected):
                    logger.warning(
                        f"Skipping row with mismatched lengths: '{original}' vs '{corrected}'"
                    )
                    continue

                total_mismatches += cls.hamming_distance(original, corrected)
                total_positions += len(original)

        if total_positions == 0:
            return None

        error_rate = total_mismatches / total_positions
        logger.debug(f"Calculated error rate: {error_rate} ({total_mismatches}/{total_positions})")
        return error_rate


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
