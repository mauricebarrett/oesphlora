import io
import os

import numpy as np
import pandas as pd
from PIL import Image
from rpy2 import robjects as r
from rpy2.robjects import StrVector, pandas2ri, r
from rpy2.robjects.vectors import ListVector
from skbio import DistanceMatrix, OrdinationResults


def plot_alpha_diversity_boxplots_with_ggplot2(
    alpha_diversity_df: pd.DataFrame,
    metadata: pd.DataFrame,
    primary_variable: str,
) -> Image.Image:
    """Plots the alpha diversity metrics for each sample, annotated with precomputed p-values,
    and saves the combined plot as a PDF.

    Args:
        alpha_diversity_df (pd.DataFrame): A DataFrame containing alpha diversity metrics for each sample.
        metadata (pd.DataFrame): A DataFrame containing the metadata for each sample.
        primary_variable (str): The primary fixed effect variable to use for plotting.

    """

    # Transfrom alpha_diversity_df to long format
    melted_df = alpha_diversity_df.melt(
        id_vars=["sample_id"],
        value_vars=[
            "Observed OTUs",
            "Simpson's dominance index",
            "Simpson's diversity index",
            "Faith's phylogenetic diversity",
        ],
        var_name="alpha_diversity_metric",
        value_name="value",
    )
    # Merge melted_df with metadata
    merged_df = melted_df.merge(
        metadata,
        how="left",
        left_on="sample_id",
        right_index=True,
    )

    # Activate pandas-to-R conversion
    pandas2ri.activate()

    # Convert pandas dataframes to R dataframes
    r_merged_df = pandas2ri.py2rpy(merged_df)

    # Assigning converted dataframes to R's
    r.assign("merged_df", r_merged_df)
    r.assign("primary_variable", primary_variable)

    # Construct the path
    # Get current path
    current_script_dir = os.path.dirname(os.path.abspath(__file__))
    # Construct the path to the R script
    r_script_path = os.path.join(current_script_dir, "../R/alpha_diversity_boxplots.R")
    # Ensure the path is absolute
    r_script_path = os.path.abspath(r_script_path)

    # Read the R code from the provided script
    with open(r_script_path, "r") as file:
        r_code = file.read()

    # Execute the R code
    temp_file_path = r(r_code)[0]

    # Read the image into Python as a PIL Image object
    with open(temp_file_path, "rb") as f:
        img = Image.open(io.BytesIO(f.read()))

    # Clean up the temporary file
    os.remove(temp_file_path)

    return img


def plot_nmds_with_ggplot2(
    nmds_coordinates: pd.DataFrame,
    metadata: pd.DataFrame,
    primary_variable: str,
    stress_value: float,
    permanova_results: dict[str, float],
) -> Image.Image:
    """
    Save NMDS coordinates and metadata, then call the R script.

    Args:
        nmds_coordinates (pd.DataFrame): DataFrame containing NMDS coordinates with samples as rows.
        metadata (pd.DataFrame): Metadata for samples.
        primary_variable (str): Column in metadata used for coloring the plot.
        stress_value (float): NMDS stress value indicating the goodness of fit.
        permanova_results (dict[str, float]): PERMANOVA results
    """

    # Activate pandas-to-R conversion
    pandas2ri.activate()

    # Convert pandas dataframes to R dataframes
    r_nmds_coordinates = pandas2ri.py2rpy(nmds_coordinates)
    r_metadata = pandas2ri.py2rpy(metadata)

    # Convert to R ListVector
    r_permanova_results = ListVector(permanova_results)

    # Assigning converted dataframes to R's
    r.assign("nmds_coordinates_df", r_nmds_coordinates)
    r.assign("metadata_df", r_metadata)
    r.assign("primary_variable", primary_variable)
    r.assign("stress_value", stress_value)
    r.assign("permanova_results", r_permanova_results)

    # Construct the path
    # Get current path
    current_script_dir = os.path.dirname(os.path.abspath(__file__))
    # Construct the path to the R script
    r_script_path = os.path.join(current_script_dir, "../R/nmds_plot.R")
    # Ensure the path is absolute
    r_script_path = os.path.abspath(r_script_path)

    # Read the R code from the provided script
    with open(r_script_path, "r") as file:
        r_code = file.read()

    # Execute the R code
    temp_file_path = r(r_code)[0]

    # Read the image into Python as a PIL Image object
    with open(temp_file_path, "rb") as f:
        img = Image.open(io.BytesIO(f.read()))

    # Clean up the temporary file
    os.remove(temp_file_path)

    return img


def add_lowest_rank_label(df):
    taxonomic_ranks = ["phylum", "class", "order", "family", "genus", "species"]

    # Determine the lowest classified rank
    df["lowest_classified"] = df.apply(
        lambda row: next(
            (
                row[rank]
                for rank in reversed(taxonomic_ranks)
                if row[rank] != "unclassified"
            ),
            "unclassified",
        ),
        axis=1,
    )

    # Create 'Label' by concatenating 'feature_id' with 'lowest_classified'
    df["label"] = df["feature_id"] + "=" + df["lowest_classified"]

    return df


def plot_biplot_with_ggplot2(
    ordination: OrdinationResults,
    metadata: pd.DataFrame,
    primary_variable: str,
    permanova_results: dict[str, float],
    taxonomy_table: pd.DataFrame,
) -> Image.Image:
    """


    Args:
        nmds_coordinates (pd.DataFrame): DataFrame containing NMDS coordinates with samples as rows.
        metadata (pd.DataFrame): Metadata for samples.
        primary_variable (str): Column in metadata used for coloring the plot.
        stress_value (float): NMDS stress value indicating the goodness of fit.
        permanova_results (dict[str, float])
    Returns:

    """

    # Merge ordination.samples with metadata using the indices
    ord_df = ordination.samples[["PC1", "PC2"]].join(
        metadata[[primary_variable]], how="inner"
    )

    # Pull out from ordination data explained variance
    explained_variance = ordination.proportion_explained
    # Convert Pandas Series to a dictionary
    explained_variance_dict = explained_variance.to_dict()

    features = ordination.features[["PC1", "PC2"]]
    norms = np.sqrt(
        (features["PC1"] ** 2 * explained_variance["PC1"])
        + (features["PC2"] ** 2 * explained_variance["PC2"])
    )
    top_features = norms.nlargest(12).index
    top_features_df = features.loc[top_features]

    # Reset the index
    top_features_df.reset_index(inplace=True)

    # Rename index column to  feature_id
    top_features_df.rename(columns={"index": "feature_id"}, inplace=True)

    # Add taxanomic information
    top_feats_taxa_df = top_features_df.merge(
        taxonomy_table, how="left", left_on="feature_id", right_index=True
    )

    # Add lowest rank label
    top_feats_taxa_df = add_lowest_rank_label(top_feats_taxa_df)

    # Activate pandas-to-R conversion
    pandas2ri.activate()

    # Convert pandas dataframes to R dataframes
    r_ord_df = pandas2ri.py2rpy(ord_df)
    r_top_feats_taxa_df = pandas2ri.py2rpy(top_feats_taxa_df)

    # Convert to R ListVector
    r_permanova_results = ListVector(permanova_results)
    r_explained_variance = ListVector(explained_variance_dict)

    # Assigning converted dataframes to R's
    r.assign("ord_df", r_ord_df)
    r.assign("top_feats_taxa_df", r_top_feats_taxa_df)
    r.assign("primary_variable", primary_variable)
    r.assign("permanova_results", r_permanova_results)
    r.assign("explained_variance", r_explained_variance)

    # Construct the path
    # Get current path
    current_script_dir = os.path.dirname(os.path.abspath(__file__))
    # Construct the path to the R script
    r_script_path = os.path.join(current_script_dir, "../R/biplot.R")
    # Ensure the path is absolute
    r_script_path = os.path.abspath(r_script_path)

    # Read the R code from the provided script
    with open(r_script_path, "r") as file:
        r_code = file.read()

    # Execute the R code
    temp_file_path = r(r_code)[0]

    # Read the image into Python as a PIL Image object
    with open(temp_file_path, "rb") as f:
        img = Image.open(io.BytesIO(f.read()))

    # Clean up the temporary file
    os.remove(temp_file_path)

    return img
