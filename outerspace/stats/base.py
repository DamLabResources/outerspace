"""Base classes for UMI statistics calculations.

This module provides abstract base classes for different types of UMI statistics
calculations including single-sample, pairwise, and differential analysis.
"""

import csv
import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, List, Optional, TypeVar, Union

from ..umi import UMI, UmiCollection

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

T = TypeVar("T")


class BaseStatistic(ABC):
    """Base class for all statistics calculations.

    This abstract base class provides a common interface for all UMI statistics
    calculations with validation, normalization, and execution methods.
    """

    def __init__(self, umi: UMI, **kwargs: Any) -> None:
        """Initialize the statistic calculator.

        Parameters
        ----------
        umi : UMI
            UMI object to calculate statistics on
        **kwargs : Any
            Additional arguments specific to the statistic
        """
        self.umi = umi
        logger.debug(f"Initialized {self.__class__.__name__} with UMI object")

    @abstractmethod
    def run(self) -> T:
        """Run the statistical calculation.

        Returns
        -------
        T
            Calculated statistic value

        Notes
        -----
        This method must be implemented by all subclasses to perform the
        actual statistical calculation.
        """
        pass

    def validate(self) -> None:
        """Validate input parameters.

        This method can be overridden by subclasses to perform validation
        of input parameters before calculation.

        Raises
        ------
        ValueError
            If validation fails
        """
        pass

    def normalize(self) -> None:
        """Normalize data before calculation if needed.

        This method can be overridden by subclasses to perform data
        normalization before statistical calculation.
        """
        pass

    @classmethod
    def calculate(cls, umi: UMI, **kwargs: Any) -> T:
        """Create an instance and run the calculation.

        Parameters
        ----------
        umi : UMI
            UMI object to calculate statistics on
        **kwargs : Any
            Additional arguments specific to the statistic

        Returns
        -------
        T
            Calculated statistic value
        """
        instance = cls(umi, **kwargs)
        instance.validate()
        instance.normalize()
        return instance.run()

    @classmethod
    def calculate_collection(
        cls, collection: UmiCollection, **kwargs: Any
    ) -> Dict[str, T]:
        """Calculate statistics for each UMI in a collection.

        Parameters
        ----------
        collection : UmiCollection
            UmiCollection object containing multiple UMIs
        **kwargs : Any
            Additional arguments specific to the statistic

        Returns
        -------
        Dict[str, T]
            Dictionary mapping sample names to calculated statistics
        """
        results = {}
        for sample_name, umi in collection.umis.items():
            logger.debug(f"Calculating {cls.__name__} for sample {sample_name}")
            results[sample_name] = cls.calculate(umi, **kwargs)
        return results

    @classmethod
    def _from_step(cls, input_file: str, sep: str = ",", **step_params: Any) -> T:
        """Create statistic from step parameters and CSV file.

        This method provides a generic interface for calculating statistics
        from configuration step parameters. Subclasses should override this
        method to implement their specific parameter handling.

        Parameters
        ----------
        input_file : str
            Path to input CSV file
        sep : str, default=','
            CSV separator
        **step_params : Any
            Step-specific parameters from configuration

        Returns
        -------
        T
            Calculated statistic value

        Raises
        ------
        NotImplementedError
            If subclass doesn't implement this method
        """
        raise NotImplementedError(
            f"{cls.__name__} must implement _from_step classmethod"
        )


class BasePairwiseStatistic(ABC):
    """Base class for pairwise statistics.

    This abstract base class provides a common interface for pairwise UMI
    statistics calculations between two UMI objects.
    """

    def __init__(self, umi1: UMI, umi2: UMI, **kwargs: Any) -> None:
        """Initialize the pairwise statistic calculator.

        Parameters
        ----------
        umi1 : UMI
            First UMI object
        umi2 : UMI
            Second UMI object
        **kwargs : Any
            Additional arguments specific to the statistic
        """
        self.umi1 = umi1
        self.umi2 = umi2
        self._result: Optional[Dict[str, T]] = None
        logger.debug(f"Initialized {self.__class__.__name__} with two UMI objects")

    @abstractmethod
    def run(self) -> T:
        """Run the pairwise statistical calculation.

        Returns
        -------
        T
            Calculated statistic value

        Notes
        -----
        This method must be implemented by all subclasses to perform the
        actual pairwise statistical calculation.
        """
        pass

    def validate(self) -> None:
        """Validate input parameters.

        This method can be overridden by subclasses to perform validation
        of input parameters before calculation.

        Raises
        ------
        ValueError
            If validation fails
        """
        pass

    def normalize(self) -> None:
        """Normalize data before calculation if needed.

        This method can be overridden by subclasses to perform data
        normalization before statistical calculation.
        """
        pass

    @classmethod
    def calculate(cls, umi1: UMI, umi2: UMI, **kwargs: Any) -> T:
        """Create an instance and run the calculation.

        Parameters
        ----------
        umi1 : UMI
            First UMI object
        umi2 : UMI
            Second UMI object
        **kwargs : Any
            Additional arguments specific to the statistic

        Returns
        -------
        T
            Calculated statistic value
        """
        instance = cls(umi1, umi2, **kwargs)
        instance.validate()
        instance.normalize()
        return instance.run()

    @classmethod
    def calculate_collection(
        cls, collection: UmiCollection, **kwargs: Any
    ) -> Dict[tuple[str, str], T]:
        """Calculate pairwise statistics for all pairs of UMIs in a collection.

        Parameters
        ----------
        collection : UmiCollection
            UmiCollection object containing multiple UMIs
        **kwargs : Any
            Additional arguments specific to the statistic

        Returns
        -------
        Dict[tuple[str, str], T]
            Dictionary mapping (sample1, sample2) tuples to calculated statistics
        """
        results = {}
        samples = list(collection.umis.keys())
        logger.info(f"Calculating pairwise statistics for {len(samples)} samples")

        for i, sample1 in enumerate(samples):
            for sample2 in samples[i + 1 :]:
                logger.debug(f"Calculating {cls.__name__} for {sample1} vs {sample2}")
                results[(sample1, sample2)] = cls.calculate(
                    collection.umis[sample1], collection.umis[sample2], **kwargs
                )
        return results


class BaseDifferentialStatistic(ABC):
    """Base class for differential analysis statistics.

    This abstract base class provides a common interface for differential
    analysis statistics across multiple UMI objects with group assignments.
    """

    def __init__(self, umis: List[UMI], groups: List[str], **kwargs: Any) -> None:
        """Initialize the differential statistic calculator.

        Parameters
        ----------
        umis : List[UMI]
            List of UMI objects
        groups : List[str]
            List of group labels corresponding to UMI objects
        **kwargs : Any
            Additional arguments specific to the statistic

        Raises
        ------
        ValueError
            If number of UMIs doesn't match number of groups
        """
        if len(umis) != len(groups):
            raise ValueError("Number of UMIs must match number of groups")
        self.umis = umis
        self.groups = groups
        self._result: Optional[Dict[str, T]] = None
        logger.debug(
            f"Initialized {self.__class__.__name__} with {len(umis)} UMIs in {len(set(groups))} groups"
        )

    @abstractmethod
    def run(self) -> T:
        """Run the differential statistical calculation.

        Returns
        -------
        T
            Calculated statistic value

        Notes
        -----
        This method must be implemented by all subclasses to perform the
        actual differential statistical calculation.
        """
        pass

    def validate(self) -> None:
        """Validate input parameters.

        This method can be overridden by subclasses to perform validation
        of input parameters before calculation.

        Raises
        ------
        ValueError
            If validation fails
        """
        pass

    def normalize(self) -> None:
        """Normalize data before calculation if needed.

        This method can be overridden by subclasses to perform data
        normalization before statistical calculation.
        """
        pass

    @classmethod
    def calculate(cls, umis: List[UMI], groups: List[str], **kwargs: Any) -> T:
        """Create an instance and run the calculation.

        Parameters
        ----------
        umis : List[UMI]
            List of UMI objects
        groups : List[str]
            List of group labels corresponding to UMI objects
        **kwargs : Any
            Additional arguments specific to the statistic

        Returns
        -------
        T
            Calculated statistic value
        """
        instance = cls(umis, groups, **kwargs)
        instance.validate()
        instance.normalize()
        return instance.run()

    @classmethod
    def calculate_collection(
        cls, collection: UmiCollection, groups: Dict[str, str], **kwargs: Any
    ) -> T:
        """Calculate differential statistics for a collection of UMIs.

        Parameters
        ----------
        collection : UmiCollection
            UmiCollection object containing multiple UMIs
        groups : Dict[str, str]
            Dictionary mapping sample names to group labels
        **kwargs : Any
            Additional arguments specific to the statistic

        Returns
        -------
        T
            Calculated statistic value

        Raises
        ------
        ValueError
            If any sample in collection is not in groups
        """
        # Validate that all samples have group assignments
        missing = set(collection.umis.keys()) - set(groups.keys())
        if missing:
            raise ValueError(f"Samples missing group assignments: {missing}")

        # Create ordered lists of UMIs and groups
        umis = []
        group_labels = []
        for sample in collection.umis.keys():
            umis.append(collection.umis[sample])
            group_labels.append(groups[sample])

        logger.info(
            f"Calculating differential statistics for {len(umis)} samples in {len(set(group_labels))} groups"
        )
        return cls.calculate(umis, group_labels, **kwargs)


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
