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
from ..nearest import NearestUMIFinder
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
    including count access.
    """

    def __init__(
        self,
        umi: UMI,
        use_corrected: bool = True,
        **kwargs: Any,
    ) -> None:
        """Initialize UMI stats calculator.

        Parameters
        ----------
        umi : UMI
            UMI object to calculate statistics on
        use_corrected : bool, default=True
            If True, use corrected counts. If False, use original counts.
        **kwargs : Any
            Additional arguments passed to parent class
        """
        super().__init__(umi, **kwargs)
        self.use_corrected = use_corrected

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
    ) -> None:
        """Initialize Gini coefficient calculator.

        Parameters
        ----------
        umi : UMI
            UMI object to calculate statistics on
        use_corrected : bool, default=True
            If True, use corrected counts. If False, use original counts.
        """
        super().__init__(umi, use_corrected)

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
        result = GiniCoefficient.calculate_gini(list(self._counts.values()))
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
            - use_corrected: bool (optional) - Whether to use corrected counts

        Returns
        -------
        Optional[float]
            Calculated Gini coefficient
        """
        key_column = step_params.get("key_column")
        barcode_column = step_params.get("barcode_column")
        use_corrected = step_params.get("use_corrected", False)  # No correction needed, data is pre-collapsed

        if not key_column:
            raise ValueError("key_column is required for GiniCoefficient")

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

        return cls.calculate(umi, use_corrected=use_corrected)


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
        """
        super().__init__(umi, use_corrected)
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
        result = ShannonDiversity.calculate_shannon(list(self._counts.values()), self.base)
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
            - use_corrected: bool (optional) - Whether to use corrected counts
            - base: float (optional) - Base for logarithm (default 2.0)

        Returns
        -------
        Optional[float]
            Calculated Shannon diversity
        """
        key_column = step_params.get("key_column")
        barcode_column = step_params.get("barcode_column")
        use_corrected = step_params.get("use_corrected", False)
        base = step_params.get("base", 2.0)

        if not key_column:
            raise ValueError("key_column is required for ShannonDiversity")

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

        return cls.calculate(umi, use_corrected=use_corrected, base=base)


class SimpsonDiversity(UMIStats):
    """Calculate Simpson's diversity index for UMI counts.

    Simpson's diversity index measures the probability that two randomly
    selected UMIs are different. Higher values indicate greater diversity.
    """

    def __init__(
        self,
        umi: UMI,
        use_corrected: bool = True,
    ) -> None:
        """Initialize Simpson's diversity calculator.

        Parameters
        ----------
        umi : UMI
            UMI object to calculate statistics on
        use_corrected : bool, default=True
            If True, use corrected counts. If False, use original counts.
        """
        super().__init__(umi, use_corrected)

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
        result = SimpsonDiversity.calculate_simpson(list(self._counts.values()))
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
            - use_corrected: bool (optional) - Whether to use corrected counts

        Returns
        -------
        Optional[float]
            Calculated Simpson diversity
        """
        key_column = step_params.get("key_column")
        barcode_column = step_params.get("barcode_column")
        use_corrected = step_params.get("use_corrected", False)

        if not key_column:
            raise ValueError("key_column is required for SimpsonDiversity")

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

        return cls.calculate(umi, use_corrected=use_corrected)


def _parse_q_parameter(q_param: Any) -> List[float]:
    """Parse q parameter into list of float values.
    
    Parameters
    ----------
    q_param : Any
        Can be a number, a string word, or a comma-separated string
        
    Returns
    -------
    List[float]
        List of q values to calculate
        
    Notes
    -----
    Supported string keywords:
    - "richness": q=0
    - "shannon": q=1
    - "simpson": q=2
    """
    # Mapping of keywords to q values
    keyword_map = {
        "richness": 0.0,
        "shannon": 1.0,
        "simpson": 2.0,
    }
    
    # If it's already a number, return it as a list
    if isinstance(q_param, (int, float)):
        return [float(q_param)]
    
    # If it's a string, parse it
    if isinstance(q_param, str):
        # Split by comma in case it's a list
        parts = [p.strip() for p in q_param.split(",")]
        q_values = []
        
        for part in parts:
            # Try to convert to float first
            try:
                q_values.append(float(part))
            except ValueError:
                # If not a number, check if it's a keyword
                part_lower = part.lower()
                if part_lower in keyword_map:
                    q_values.append(keyword_map[part_lower])
                else:
                    raise ValueError(
                        f"Invalid q parameter '{part}'. "
                        f"Must be a number or one of: {', '.join(keyword_map.keys())}"
                    )
        
        return q_values
    
    raise ValueError(f"Invalid q parameter type: {type(q_param)}")


class HillNumber(UMIStats):
    """Calculate Hill numbers for diversity analysis.
    
    Hill numbers are a mathematically unified family of diversity indices
    differing among themselves only by an exponent q that determines their
    sensitivity to species relative abundances.
    """

    def __init__(
        self,
        umi: UMI,
        q: float = 1.0,
        use_corrected: bool = True,
    ) -> None:
        """Initialize Hill number calculator.

        Parameters
        ----------
        umi : UMI
            UMI object to calculate statistics on
        q : float, default=1.0
            Order of the Hill number. Common values:
            - q=0: Species richness (count of species)
            - q=1: Exponential of Shannon entropy
            - q=2: Inverse Simpson concentration
            - q→∞: Inverse of max proportional abundance
        use_corrected : bool, default=True
            If True, use corrected counts. If False, use original counts.
        """
        super().__init__(umi, use_corrected)
        self.q = q

    @staticmethod
    def calculate_hill(counts: Sequence[float], q: float) -> Optional[float]:
        """Calculate Hill number for given q parameter.

        Parameters
        ----------
        counts : Sequence[float]
            Sequence of count values
        q : float
            Order of the Hill number

        Returns
        -------
        Optional[float]
            Hill number or None if calculation is not possible

        Notes
        -----
        Hill numbers are a generalized diversity measure calculated as:
        
        For q ≠ 1:
            D_q = (sum(p_i^q))^(1/(1-q))
            
        For q = 1 (limit as q→1):
            D_1 = exp(-sum(p_i * ln(p_i))) = exp(H)
            where H is Shannon entropy
            
        Special cases:
        - q=0: Species richness (count of non-zero elements)
        - q=1: Exponential of Shannon entropy
        - q=2: Inverse Simpson concentration (1 / sum(p_i^2))
        - q→∞: Inverse of max proportional abundance (1 / max(p_i))
        
        References
        ----------
        Hill, M. O. (1973). Diversity and evenness: a unifying notation and its
        consequences. Ecology, 54(2), 427-432.
        """
        if not counts:
            return None

        # Filter out zero counts
        counts_array = np.array([c for c in counts if c > 0])
        
        if len(counts_array) == 0:
            return None
            
        total = np.sum(counts_array)
        
        if total == 0:
            return None

        # Calculate proportions
        proportions = counts_array / total

        # Special case: q = 0 (species richness)
        if q == 0:
            return float(len(counts_array))
        
        # Special case: q = 1 (exponential Shannon entropy)
        # Use limit as q→1: exp(-sum(p_i * ln(p_i)))
        elif abs(q - 1.0) < 1e-9:  # Close to 1
            shannon = -np.sum(proportions * np.log(proportions))
            return float(np.exp(shannon))
        
        # General case: q ≠ 1
        else:
            sum_p_q = np.sum(proportions ** q)
            hill = sum_p_q ** (1 / (1 - q))
            return float(hill)

    def run(self) -> Optional[float]:
        """Calculate Hill number for the UMI object.

        Returns
        -------
        Optional[float]
            Hill number or None if calculation is not possible
        """
        result = HillNumber.calculate_hill(list(self._counts.values()), self.q)
        logger.debug(f"Calculated Hill number (q={self.q}): {result}")
        return result

    @classmethod
    def _from_step(cls, input_file: str, sep: str = ",", **step_params: Any) -> Any:
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
            - use_corrected: bool (optional) - Whether to use corrected counts
            - q: float, str, or comma-separated str - Order parameter(s) (required)
              Can be:
              - A number: 0, 1, 2, 1.5, etc.
              - A keyword: "richness" (0), "shannon" (1), "simpson" (2)
              - Comma-separated: "richness, shannon, simpson, 1.5"

        Returns
        -------
        float or Dict[str, float]
            If single q value: returns float
            If multiple q values: returns dict mapping q names to values
        """
        key_column = step_params.get("key_column")
        barcode_column = step_params.get("barcode_column")
        use_corrected = step_params.get("use_corrected", False)
        q_param = step_params.get("q")

        if not key_column:
            raise ValueError("key_column is required for HillNumber")
        if q_param is None:
            raise ValueError("q parameter is required for HillNumber")

        # Parse q parameter
        q_values = _parse_q_parameter(q_param)

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

        # If single q value, return single result
        if len(q_values) == 1:
            return cls.calculate(umi, q=q_values[0], use_corrected=use_corrected)
        
        # If multiple q values, return dict
        results = {}
        for q in q_values:
            result = cls.calculate(umi, q=q, use_corrected=use_corrected)
            results[f"q={q}"] = result
        
        return results



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
        super().__init__(umi, use_corrected)
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
        super().__init__(umi, use_corrected)
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


class ErrorRate(BaseStatistic):
    """Calculate error rate by comparing original and corrected columns using alignment scoring.

    This class compares two columns from a CSV file to calculate the
    error rate based on alignment scores from NearestUMIFinder. The error rate
    is calculated as the normalized difference between the maximum possible score
    and the actual alignment score.
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
            - original_column: str - Column with original values (required)
            - corrected_column: str - Column with corrected values (required)
            - mismatch_penalty: int - Penalty for mismatches (default: -1)
            - gap_penalty: int - Penalty for gaps (default: -3)
            - match_score: int - Score for matches (default: 1)

        Returns
        -------
        Optional[float]
            Calculated error rate (normalized between 0 and 1, where 0 = perfect match)

        Notes
        -----
        The error rate is calculated as:
        error_rate = (max_possible_score - actual_score) / max_possible_score
        
        This uses the alignment scoring system from NearestUMIFinder, which accounts
        for mismatches and gaps more accurately than simple Hamming distance.
        """
        original_column = step_params.get("original_column")
        corrected_column = step_params.get("corrected_column")
        mismatch_penalty = step_params.get("mismatch_penalty", -1)
        gap_penalty = step_params.get("gap_penalty", -3)
        match_score = step_params.get("match_score", 1)

        if not original_column or not corrected_column:
            raise ValueError("Both original_column and corrected_column are required for ErrorRate")

        total_max_score = 0
        total_actual_score = 0
        num_comparisons = 0

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

                # Calculate maximum possible score (perfect match)
                max_score = len(original) * match_score
                
                # Calculate actual alignment score
                actual_score = NearestUMIFinder._calculate_alignment_score(
                    original,
                    corrected,
                    mismatch_penalty,
                    gap_penalty,
                    match_score,
                    None  # No cutoff
                )
                
                total_max_score += max_score
                total_actual_score += actual_score
                num_comparisons += 1

        if total_max_score == 0 or num_comparisons == 0:
            return None

        # Calculate normalized error rate
        error_rate = (total_max_score - total_actual_score) / total_max_score
        logger.debug(
            f"Calculated error rate: {error_rate:.6f} "
            f"(max_score={total_max_score}, actual_score={total_actual_score}, "
            f"n={num_comparisons})"
        )
        return error_rate


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
