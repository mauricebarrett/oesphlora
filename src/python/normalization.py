import numpy as np
import pandas as pd


def rclr_transform(
    counts: pd.DataFrame, qiime_orientation: bool = False
) -> pd.DataFrame:
    """
    Robust Center Log-Ratio (rCLR) transformation for sparse microbiome count data.

    Performs log transformation and centers each sample using the median of non-zero
    logged values. More robust than standard CLR for compositional data with many zeros.
    Automatically removes samples with all-zero counts.

    Expected Data Structure:
        Default (qiime_orientation=False):
            - Rows: Samples
            - Columns: Features (ASVs/OTUs/taxa)
        QIIME2 format (qiime_orientation=True):
            - Rows: Features (ASVs/OTUs/taxa)
            - Columns: Samples
            - Will be transposed internally to samples × features
        Values: Non-negative counts

    Transformation:
        rCLR[i,j] = log(counts[i,j]) - median(log(counts[i,k])) for k where counts[i,k] > 0
        Zeros remain as NaN in output.

    Args:
        counts (pd.DataFrame): Count table with non-negative values.
        qiime_orientation (bool): If True, transpose input from QIIME2 format
            (features × samples) to standard format (samples × features). Default False.

    Returns:
        pd.DataFrame: rCLR-transformed data with all-zero samples removed.
            Always returns in samples × features orientation.
            Zeros become NaN, non-zeros are log-transformed and sample-centered.

    Raises:
        ValueError: If input contains negative values.

    Example:
        >>> # Standard orientation (samples × features)
        >>> rclr_data = rclr_transform(asv_counts)
        >>>
        >>> # QIIME2 orientation (features × samples)
        >>> rclr_data = rclr_transform(qiime_table, qiime_orientation=True)
    """
    if (counts.values < 0).any():
        raise ValueError("Counts must be non-negative.")

    # Handle QIIME2 orientation by transposing
    if qiime_orientation:
        counts = counts.T

    # log of nonzeros; zeros -> NaN
    counts_with_nan = counts.replace(0, np.nan)
    log_nonzero = pd.DataFrame(
        np.log(counts_with_nan),
        index=counts.index,
        columns=counts.columns,
    )

    # robust center per sample (row-wise median of logged nonzeros)
    row_center = np.nanmedian(log_nonzero.values, axis=1)

    # subtract the center from each row
    rclr = log_nonzero.sub(row_center, axis=0)

    # always drop all-zero samples (they have NaN row_center)
    rclr = rclr.loc[~np.isnan(row_center)]

    return rclr


def presence_absence_transform(
    counts: pd.DataFrame, qiime_orientation: bool = False
) -> pd.DataFrame:
    """
    Presence-absence transformation for microbiome count data.

    Converts counts to binary presence (1) or absence (0) based on non-zero values.
    Automatically removes samples with all-zero counts.

    Args:
        counts (pd.DataFrame): Count table with non-negative values.
        qiime_orientation (bool): If True, transpose input from QIIME2 format
            (features × samples) to standard format (samples × features). Default False.

    Returns:
        pd.DataFrame: Presence-absence transformed data with all-zero samples removed.
            Always returns in samples × features orientation.
    """
    if qiime_orientation:
        counts = counts.T

    # Convert to presence-absence
    pa_data = (counts > 0).astype(int)

    # Drop all-zero samples
    pa_data = pa_data.loc[pa_data.sum(axis=1) > 0]

    return pa_data
