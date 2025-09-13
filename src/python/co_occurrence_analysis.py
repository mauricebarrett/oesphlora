import os

import pandas as pd
from rpy2.robjects import pandas2ri, r  # type: ignore
from rpy2.robjects.vectors import ListVector  # type: ignore


def prevalance_filtering(
    feature_table_df: pd.DataFrame, prevalence_threshold: float
) -> pd.DataFrame:
    """
    Remove features from feature table with prevalence below threshold.

    Filters out taxonomic features (ASVs/OTUs/Taxa) that are present in fewer
    samples than the specified percentage threshold. This removes rare
    features that may represent noise or sequencing artifacts.

    Args:
        feature_table_df (pd.DataFrame): Feature abundance table with features
            as rows and samples as columns. Values should be non-negative.
        prevalence_threshold (float): Minimum prevalence percentage (0-100).
            Features must be present (>0) in at least this percentage of samples.

    Returns:
        pd.DataFrame: Filtered table containing only features meeting the
            prevalence threshold.

    Example:
        >>> # Keep features present in ≥10% of samples
        >>> filtered_table = prevalance_filtering(asv_table, 10.0)
        >>> # Original: 1000 features → Filtered: 450 features

    Note:
        Prevalence = (number of samples where feature > 0) / total samples * 100
    """
    feature_prevalence = (
        (feature_table_df > 0).sum(axis=1) / feature_table_df.shape[1] * 100
    )

    # print the number of features before filtering
    print(f"Number of features before filtering: {feature_table_df.shape[0]}")

    # Filter the dataframe FIRST, then print the shape
    filtered_df = feature_table_df.loc[feature_prevalence >= prevalence_threshold, :]

    # print the number of features after filtering
    print(f"Number of features after filtering: {filtered_df.shape[0]}")
    print(f"Shape of feature_table_df after filtering: {filtered_df.shape}")

    # Return the filtered dataframe
    return filtered_df


def perform_co_occurrence_network(
    taxa_table_df: pd.DataFrame,
    title_dict: dict,
    output_dir: str,
) -> None:
    """
    Create a co-occurrence network from the OTU table and metadata.
    """
    # create output directory if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # filter taxa table by prevalence
    taxa_table_df = prevalance_filtering(taxa_table_df, 20)

    # Activate pandas-to-R conversion
    pandas2ri.activate()

    # Convert pandas dataframes to R dataframes
    r_taxa_table_df = pandas2ri.py2rpy(taxa_table_df)
    r_title_dict = ListVector(title_dict)

    # Assign variables to R global environment (simpler approach)
    r.assign("taxa_table_df", r_taxa_table_df)
    r.assign("title_list", r_title_dict)
    r.assign("output_dir", output_dir)

    # Construct the path to R script
    current_script_dir = os.path.dirname(os.path.abspath(__file__))
    r_script_path = os.path.join(current_script_dir, "../R/co_occurrence_analysis.R")
    r_script_path = os.path.abspath(r_script_path)

    # Read and execute the R code
    with open(r_script_path, "r") as file:
        r_code = file.read()

    # Execute the R code
    r(r_code)
