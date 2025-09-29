"""UMI clustering and correction functionality.

This module provides classes and functions for handling Unique Molecular Identifiers (UMIs),
including clustering, correction, and management across multiple samples. It supports
various clustering methods and can work with CSV files and pandas DataFrames.
"""

import logging
from typing import Dict, List, Optional, Union, Tuple
from collections import Counter
from pathlib import Path
import csv
import sys

from umi_tools import UMIClusterer

# Set up logging
logger = logging.getLogger(__name__)

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"


def _cluster_umis(
    counts: Dict[bytes, int],
    mismatches: int = 2,
    method: str = "adjacency",
    allowed_list: Optional[List[bytes]] = None,
) -> Tuple[Dict[bytes, int], Dict[bytes, bytes]]:
    """Cluster UMIs and return collapsed counts and mapping.

    This function performs UMI clustering using the specified method and returns
    both the collapsed count dictionary and the mapping from original to
    representative UMIs.

    Parameters
    ----------
    counts : Dict[bytes, int]
        Dictionary mapping UMIs to their counts
    mismatches : int, default=2
        Number of mismatches allowed for clustering
    method : str, default='adjacency'
        Clustering method to use ('adjacency', 'cluster', or 'directional')
    allowed_list : Optional[List[bytes]], default=None
        Optional list of allowed UMIs. If provided, highest priority UMI in cluster
        will be used as representative

    Returns
    -------
    Tuple[Dict[bytes, int], Dict[bytes, bytes]]
        Tuple containing:
        - collapsed_counts: Dictionary mapping representative UMIs to total counts
        - mapping: Dictionary mapping original UMIs to representative UMIs

    Notes
    -----
    If no mismatches are allowed (mismatches=0), each UMI maps to itself.
    Clusters are processed in order of size (largest first) to ensure
    consistent representative selection.
    """
    if not counts:
        return {}, {}

    # If no mismatches allowed, each barcode maps to itself
    if mismatches == 0:
        return counts.copy(), {bc: bc for bc in counts.keys()}

    # Each UMI should only occur once
    allowed_set = set(allowed_list) if allowed_list else set()

    # Create clusters
    clusterer = UMIClusterer(cluster_method=method)
    clusters = clusterer(counts, mismatches)

    # Sort clusters by size (largest first for consistent processing)
    clusters = sorted(clusters, key=lambda x: _cluster_size(counts, x), reverse=True)

    # Create mapping
    mapping = {}
    collapsed_counts = {}

    for cluster in clusters:
        if allowed_list:
            # Find highest priority UMI in cluster that's in allowed_list
            key = None
            # Find the UMI in the cluster that is in the allowed list
            # and has the highest count
            for bc in cluster:
                if bc in allowed_set:
                    if key is None or counts[bc] > counts[key]:
                        key = bc
            if key is None:
                # If no allowed UMIs in cluster, use most common
                key = max(cluster, key=lambda x: counts[x])
        else:
            # Use most common UMI as representative
            key = max(cluster, key=lambda x: counts[x])

        # Map all UMIs in cluster to representative
        for item in cluster:
            mapping[item] = key
            collapsed_counts[key] = collapsed_counts.get(key, 0) + counts[item]

        # Remove the key from the allowed set
        allowed_set.discard(key)

    return collapsed_counts, mapping


def _cluster_size(counts: Dict[bytes, int], cluster: List[bytes]) -> int:
    """Calculate the total size of a cluster.

    Parameters
    ----------
    counts : Dict[bytes, int]
        Dictionary mapping UMIs to their counts
    cluster : List[bytes]
        List of UMIs in the cluster

    Returns
    -------
    int
        Total count of all UMIs in the cluster
    """
    return sum(counts[bc] for bc in cluster)



class UMI:
    """Class for handling UMI clustering and correction.

    This class provides functionality for UMI clustering, correction, and
    management. It supports various clustering methods and can work with
    both individual UMIs and bulk data from CSV files.
    """

    def __init__(
        self, mismatches: int = 2, method: str = "adjacency", correct: bool = True
    ) -> None:
        """Initialize UMI clusterer.

        Parameters
        ----------
        mismatches : int, default=2
            Number of mismatches allowed for clustering
        method : str, default='adjacency'
            Clustering method to use ('cluster', 'adjacency', or 'directional')
        correct : bool, default=True
            Whether to perform clustering (default: True)
        """
        self.clusterer = UMIClusterer(cluster_method=method)
        self.mismatches = mismatches
        self.method = method  # Store the method name
        self._counts: Dict[bytes, int] = {}
        self._mapping: Dict[bytes, bytes] = {}
        self._corrected_counts: Optional[Dict[bytes, int]] = None
        self.correct = correct

        logger.debug(
            f"Initialized UMI clusterer: {method} method, {mismatches} mismatches"
        )

    def consume(self, umi: Union[str, bytes], n: int = 1) -> None:
        """Add a UMI to the counts dictionary.

        Parameters
        ----------
        umi : Union[str, bytes]
            UMI sequence to add
        n : int, default=1
            Number of times to add the UMI
        """
        if isinstance(umi, str):
            umi = umi.encode("ascii")
        self._counts[umi] = self._counts.get(umi, 0) + n
        self._corrected_counts = None  # Reset corrected counts

    def create_mapping(
        self, allowed_list: Optional[List[Union[str, bytes]]] = None
    ) -> None:
        """Create mapping between original and corrected barcodes.

        Parameters
        ----------
        allowed_list : Optional[List[Union[str, bytes]]], default=None
            Optional list of allowed UMIs. If provided, highest priority UMI in cluster
            will be used as representative

        Notes
        -----
        If correction is disabled, each barcode maps to itself.
        """
        if not self._counts:
            logger.debug("No counts available for mapping creation")
            return

        # Convert allowed_list to bytes if provided
        if allowed_list:
            allowed_list = [
                x.encode("ascii") if isinstance(x, str) else x for x in allowed_list
            ]

        # If correction disabled, each barcode maps to itself
        if not self.correct:
            self._mapping = {bc: bc for bc in self._counts.keys()}
            self._corrected_counts = self._counts.copy()
            logger.debug("Correction disabled - using original counts")
            return

        # Use private clustering function
        self._corrected_counts, self._mapping = _cluster_umis(
            self._counts,
            self.mismatches,
            self.method,  # Use stored method name
            allowed_list,
        )

        logger.debug(f"Created mapping for {len(self._counts)} UMIs")

    def __getitem__(self, umi: Union[str, bytes]) -> bytes:
        """Get corrected UMI for a given UMI.

        Parameters
        ----------
        umi : Union[str, bytes]
            Original UMI sequence

        Returns
        -------
        bytes
            Corrected UMI sequence
        """
        if isinstance(umi, str):
            umi = umi.encode("ascii")
        return self._mapping.get(umi, umi)

    def get(
        self, umi: Union[str, bytes], default: Optional[bytes] = None
    ) -> Optional[bytes]:
        """Get corrected UMI for a given UMI with optional default.

        Parameters
        ----------
        umi : Union[str, bytes]
            Original UMI sequence
        default : Optional[bytes], default=None
            Default value if UMI not found

        Returns
        -------
        Optional[bytes]
            Corrected UMI sequence or default value
        """
        if isinstance(umi, str):
            umi = umi.encode("ascii")
        return self._mapping.get(umi, default)

    @property
    def corrected_counts(self) -> Dict[bytes, int]:
        """Get counts of corrected barcodes.

        Returns
        -------
        Dict[bytes, int]
            Dictionary mapping corrected UMIs to their counts

        Notes
        -----
        If no mapping exists, returns original counts. Otherwise, creates
        corrected counts from the mapping.
        """
        if self._corrected_counts is None:
            # If no mapping exists, use original counts
            if not self._mapping:
                self._corrected_counts = self._counts.copy()
            else:
                # Create corrected counts from mapping
                self._corrected_counts = {}
                for orig_bc, count in self._counts.items():
                    corrected = self._mapping.get(orig_bc, orig_bc)
                    self._corrected_counts[corrected] = (
                        self._corrected_counts.get(corrected, 0) + count
                    )
        return self._corrected_counts

    @classmethod
    def from_csv(
        cls,
        filepath: Union[str, Path],
        column: str,
        mismatches: int = 2,
        method: str = "adjacency",
        sep: str = ",",
        correct: bool = True,
        count_column: Optional[str] = None,
        scale: Optional[float] = None,
    ) -> "UMI":
        """Create UMI object from CSV file.

        Parameters
        ----------
        filepath : Union[str, Path]
            Path to CSV file
        column : str
            Column containing UMIs
        mismatches : int, default=2
            Number of mismatches allowed for clustering
        method : str, default='adjacency'
            Clustering method to use
        sep : str, default=','
            CSV separator
        correct : bool, default=True
            If True, perform clustering. If False, assume data is pre-clustered.
        count_column : Optional[str], default=None
            Optional column containing pre-counted values
        scale : Optional[float], default=None
            Optional scale factor for normalized values

        Returns
        -------
        UMI
            UMI object with loaded data

        Raises
        ------
        ValueError
            If required columns are not found in the CSV file
        """
        umi = cls(mismatches=mismatches, method=method, correct=correct)

        logger.info(f"Loading UMI data from {filepath}")

        with open(filepath, "r") as f:
            reader = csv.DictReader(f, delimiter=sep)
            if column not in reader.fieldnames:
                raise ValueError(f"Column {column} not found in CSV file")
            if count_column and count_column not in reader.fieldnames:
                raise ValueError(f"Count column {count_column} not found in CSV file")

            for row in reader:
                value = row[column]
                if not value:  # Skip empty values
                    continue

                if count_column:
                    try:
                        count = float(row[count_column])
                        if scale:
                            count = int(round(count * scale))
                        else:
                            count = int(round(count))
                    except (ValueError, TypeError):
                        logger.warning(
                            f"Invalid count value '{row[count_column]}' for {value}"
                        )
                        continue

                    # Add the value count times
                    umi.consume(value, count)
                else:
                    umi.consume(value)

        if correct:
            umi.create_mapping()
        else:
            # For pre-clustered data, use each UMI as its own cluster
            umi._mapping = {bc: bc for bc in umi._counts.keys()}
            umi._corrected_counts = umi._counts.copy()

        logger.info(f"Loaded {len(umi._counts)} UMIs from {filepath}")
        return umi


class UmiCollection:
    """Class for managing collections of UMI objects across samples.

    This class provides functionality for managing multiple UMI objects,
    one per sample, and performing operations across all samples such as
    collapsing UMIs and converting to various output formats.
    """

    def __init__(self, umis: Optional[Dict[str, UMI]] = None) -> None:
        """Initialize UmiCollection.

        Parameters
        ----------
        umis : Optional[Dict[str, UMI]], default=None
            Dictionary mapping sample names to UMI objects
        """
        self.umis = umis or {}
        self.logger = logging.getLogger(__name__)

    def collapse_umis(
        self,
        mismatches: int = 2,
        method: str = "adjacency",
        allowed_list: Optional[List[Union[str, bytes]]] = None,
    ) -> "UmiCollection":
        """Collapse UMIs across all samples using clustering.

        Parameters
        ----------
        mismatches : int, default=2
            Number of mismatches allowed for clustering
        method : str, default='adjacency'
            Clustering method to use
        allowed_list : Optional[List[Union[str, bytes]]], default=None
            Optional list of allowed UMIs. If provided, highest priority UMI in cluster
            will be used as representative

        Returns
        -------
        UmiCollection
            New UmiCollection with collapsed UMIs

        Notes
        -----
        This method combines all UMIs from all samples, performs clustering,
        and creates a new collection with the collapsed results.
        """
        # Combine all UMIs from all samples
        combined_counts = {}
        for umi_obj in self.umis.values():
            for bc, count in umi_obj._counts.items():
                combined_counts[bc] = combined_counts.get(bc, 0) + count

        # Convert allowed_list to bytes if provided
        if allowed_list:
            allowed_list = [
                x.encode("ascii") if isinstance(x, str) else x for x in allowed_list
            ]

        # Use private clustering function
        collapsed_counts, mapping = _cluster_umis(
            combined_counts, mismatches, method, allowed_list
        )

        # Create new collection with collapsed counts
        new_umis = {}
        for sample_name, umi_obj in self.umis.items():
            new_umi = UMI(mismatches=mismatches, method=method, correct=False)
            new_umi._mapping = mapping
            new_umi._counts = umi_obj._counts.copy()

            # Update corrected counts using the mapping
            new_umi._corrected_counts = {}
            for orig_bc, count in umi_obj._counts.items():
                corrected = mapping.get(orig_bc, orig_bc)
                new_umi._corrected_counts[corrected] = (
                    new_umi._corrected_counts.get(corrected, 0) + count
                )

            new_umis[sample_name] = new_umi

        logger.info(f"Collapsed UMIs across {len(self.umis)} samples")
        return UmiCollection(new_umis)

    @classmethod
    def from_csvs(
        cls,
        filepaths: List[Union[str, Path]],
        column: str,
        sample_names: Optional[List[str]] = None,
        mismatches: int = 0,
        method: str = "adjacency",
        sep: str = ",",
        count_column: Optional[str] = None,
    ) -> "UmiCollection":
        """Create UmiCollection from multiple CSV files.

        Parameters
        ----------
        filepaths : List[Union[str, Path]]
            List of paths to CSV files
        column : str
            Column containing UMIs
        sample_names : Optional[List[str]], default=None
            Optional list of sample names. If not provided, uses basenames
        mismatches : int, default=0
            Number of mismatches allowed for clustering
        method : str, default='adjacency'
            Clustering method to use
        sep : str, default=','
            CSV separator
        count_column : Optional[str], default=None
            Optional column containing counts

        Returns
        -------
        UmiCollection
            UmiCollection object with loaded data

        Raises
        ------
        ValueError
            If number of sample names doesn't match number of files, or if
            required columns are not found in CSV files
        """

        if sample_names and len(sample_names) != len(filepaths):
            raise ValueError("Number of sample names must match number of files")

        umis = {}
        for i, filepath in enumerate(filepaths):
            # Use provided sample name or basename
            sample_name = sample_names[i] if sample_names else Path(filepath).stem
            logger.info(
                f"Processing file {i+1}/{len(filepaths)}: {filepath} as sample {sample_name}"
            )

            # Create UMI object for this file
            umi = UMI(mismatches=mismatches, method=method, correct=False)

            # Read CSV and add counts
            with open(filepath, "r") as f:
                reader = csv.DictReader(f, delimiter=sep)
                if column not in reader.fieldnames:
                    raise ValueError(f"Column {column} not found in {filepath}")
                if count_column and count_column not in reader.fieldnames:
                    raise ValueError(
                        f"Count column {count_column} not found in {filepath}"
                    )

                for row in reader:
                    value = row[column]
                    if not value:  # Skip empty values
                        continue

                    if count_column:
                        try:
                            count = int(float(row[count_column]))
                            umi.consume(value, count)
                        except (ValueError, TypeError):
                            logger.warning(
                                f"Invalid count value '{row[count_column]}' for {value}"
                            )
                            continue
                    else:
                        umi.consume(value)

            umis[sample_name] = umi

        logger.info(f"Loaded UmiCollection from {len(filepaths)} files")
        return cls(umis)

    @classmethod
    def from_df(
        cls,
        df: "pd.DataFrame",
        sample_col: str,
        umi_col: str,
        count_col: Optional[str] = None,
    ) -> "UmiCollection":
        """Create UmiCollection from pandas DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame containing UMI data
        sample_col : str
            Column containing sample names
        umi_col : str
            Column containing UMIs
        count_col : Optional[str], default=None
            Optional column containing counts

        Returns
        -------
        UmiCollection
            UmiCollection object with loaded data
        """
        umis = {}

        for sample in df[sample_col].unique():
            logger.info(f"Processing sample: {sample}")
            sample_df = df[df[sample_col] == sample]
            umi = UMI(mismatches=0, correct=False)

            if count_col:
                for _, row in sample_df.iterrows():
                    umi.consume(row[umi_col], int(row[count_col]))
            else:
                for umi_seq in sample_df[umi_col]:
                    umi.consume(umi_seq)

            umis[sample] = umi

        logger.info(f"Loaded UmiCollection from DataFrame with {len(umis)} samples")
        return cls(umis)

    def to_df(self, format: str = "wide") -> "pd.DataFrame":
        """Convert UmiCollection to pandas DataFrame.

        Parameters
        ----------
        format : str, default='wide'
            Output format, either 'wide' or 'long'
            - wide: Each sample is a column, rows are UMIs
            - long: Three columns: sample, umi, count

        Returns
        -------
        pd.DataFrame
            DataFrame in specified format

        Raises
        ------
        ValueError
            If format is not 'wide' or 'long'
        """
        import pandas as pd

        if format not in ["wide", "long"]:
            raise ValueError("Format must be either 'wide' or 'long'")

        rows = []
        for sample, umi in self.umis.items():
            for bc, count in umi.corrected_counts.items():
                rows.append(
                    {"sample": sample, "umi": bc.decode("ascii"), "count": count}
                )
        long = pd.DataFrame(rows)

        if format == "wide":
            return pd.pivot_table(
                long, values="count", index="umi", columns="sample", fill_value=0
            )

        return long

    def write(
        self, path: Union[str, Path], sep: str = ",", format: str = "wide"
    ) -> None:
        """Write UmiCollection to CSV file.

        Parameters
        ----------
        path : Union[str, Path]
            Path to output CSV file
        sep : str, default=','
            CSV separator
        format : str, default='wide'
            Output format, either 'wide' or 'long'
        """
        logger.info(f"Writing data to {path} in {format} format")
        df = self.to_df(format=format)
        df.to_csv(path, sep=sep, index=(format == "wide"))
        logger.info(f"Successfully wrote {len(df)} rows to {path}")


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.
