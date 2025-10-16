"""Statistics module for UMI analysis"""

__copyright__ = "Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved."
__author__ = "WND"

from .base import BaseStatistic, BasePairwiseStatistic, BaseDifferentialStatistic

# Import statistics as they are implemented
from .single import (
    GiniCoefficient,
    ShannonDiversity,
    SimpsonDiversity,
    HillNumber,
    UMIRecoveryRate,
    UMIEfficiencyRate,
    ErrorRate,
)

from .multi import (
    JaccardSimilarity,
    BrayCurtisDissimilarity,
    SpearmanCorrelation,
    FoldChange,
)

from .differential import (
    SingleSampleDifferentialAbundance,
    MannWhitneyDifferentialAbundance,
    PairedTTestDifferentialAbundance,
)


__all__ = [
    "BaseStatistic",
    "BasePairwiseStatistic",
    "BaseDifferentialStatistic",
    "GiniCoefficient",
    "ShannonDiversity",
    "SimpsonDiversity",
    "HillNumber",
    "UMIRecoveryRate",
    "UMIEfficiencyRate",
    "ErrorRate",
    "JaccardSimilarity",
    "BrayCurtisDissimilarity",
    "SpearmanCorrelation",
    "FoldChange",
    "SingleSampleDifferentialAbundance",
    "MannWhitneyDifferentialAbundance",
    "PairedTTestDifferentialAbundance"
    # Add other classes as they are implemented
]


# Copyright (C) 2025, SC Barrera, R Berman, Drs DVK & WND. All Rights Reserved.