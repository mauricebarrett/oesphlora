import itertools
import os
import warnings

import numpy as np
import pandas as pd
from rpy2 import robjects
from rpy2.robjects import pandas2ri, r


def prevalance_filtering(
    otu_table: pd.DataFrame, prevalence_threshold: float
) -> pd.DataFrame:
    """Remove features from feature table with prevalence below threshold."""
    feature_prevalence = (otu_table > 0).sum(axis=1) / otu_table.shape[1] * 100

    # Remove features with prevalence below threshold
    return otu_table.loc[feature_prevalence >= prevalence_threshold, :]


def linda_daa_asv(
    asv_table: pd.DataFrame,
    metadata: pd.DataFrame,
    primary_variable: str,
    taxonomy_table: pd.DataFrame,
    output_file: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the LINDA analysis.

    Args:
        asv_table (pd.DataFrame): The OTU/feature count table.
        metadata (pd.DataFrame): The metadata.
        primary_variable (str): The column name for the variable of interest.
        taxonomy_table (pd.DataFrame): Table with taxonomic information.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: Full results and filtered significant results.
    """

    # Define path to filtered DAA results by replacing .csv extension with _filtered.csv
    filtered_output_file = output_file.replace(".csv", "_filtered.csv")

    # Skip if the filtered results file already exists
    if os.path.exists(filtered_output_file):
        print(
            f"Filtered results file {filtered_output_file} already exists. Skipping analysis."
        )
        return pd.read_csv(output_file), pd.read_csv(filtered_output_file)

    # Drop rows with missing values in the primary variable
    metadata = metadata.dropna(subset=[primary_variable])

    unique_groups = metadata[primary_variable].unique()
    pairwise_combinations = list(itertools.combinations(unique_groups, 2))

    # Create a list to store the results of the analysis
    results = []

    for group_a, group_b in pairwise_combinations:
        print(f"Processing pair: {group_a} vs {group_b}")
        # subset metadata to only include the two groups - CREATE EXPLICIT COPY
        metadata_subset = metadata[
            metadata[primary_variable].isin([group_a, group_b])
        ].copy()

        # Check if the subset has at least one sample for each group
        if metadata_subset[primary_variable].nunique() != 2:
            print(
                f"Skipping {group_a} vs {group_b}: insufficient data in metadata_subset"
            )
            continue

        # Subset OTU table base on metadata
        taxon_table_fil = asv_table[metadata_subset.index]

        taxon_table_fil.isna().any().any()

        # Skip if there are fewer than 4 samples (columns) - LINDA's correlation calculations need sufficient samples
        if taxon_table_fil.shape[1] < 4:
            print(
                f"Skipping {group_a} vs {group_b}: only {taxon_table_fil.shape[1]} samples found (minimum 4 required)"
            )
            continue

        # Suppress FutureWarning from rpy2
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning)
            pandas2ri.activate()

            # Convert pandas dataframes to R dataframes
            r_species_matrix = pandas2ri.py2rpy(taxon_table_fil)
            r_metadata = pandas2ri.py2rpy(metadata_subset)

        # Assigning converted dataframes to R's global environment
        r.assign("r_species_matrix", r_species_matrix)
        r.assign("r_metadata", r_metadata)
        formula = f"~ {primary_variable}"
        r.assign("formula", formula)

        # Retrieve the factor levels of the primary variable in R
        levels_r = r(f"levels(as.factor(r_metadata${primary_variable}))")
        levels = list(levels_r)

        r_group_a, r_group_b = levels

        # Call the LINDA function
        r(
            """
            library(MicrobiomeStat)
    

            linda_output = MicrobiomeStat::linda(
                r_species_matrix,
                r_metadata,
                formula = formula,
                feature.dat.type = 'count',
                prev.filter = 0.10,
                mean.abund.filter = 0,
                max.abund.filter = 0,
                is.winsor = TRUE,
                adaptive = TRUE,
                pseudo.cnt = 0.5,
                corr.cut = 0.1,
                n.cores = 1,
                verbose = TRUE
            )

            linda_results <- as.data.frame(linda_output$output[[1]])
            """
        )

        # Convert results to pandas dataframe - suppress warnings again
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning)
            linda_results = r["linda_results"]
            linda_results_df = pandas2ri.rpy2py(linda_results)

        # Add negative log10 p-value column
        linda_results_df["neg_log10_pvalue"] = -np.log10(linda_results_df["padj"])

        linda_results_df["direction"] = linda_results_df["log2FoldChange"].apply(
            lambda log2fc: f"Up in {r_group_a}" if log2fc > 0 else f"Up in {r_group_b}"
        )

        linda_results_df["comparison"] = f"{r_group_a} vs {r_group_b}"
        results.append(linda_results_df)

    # Combine all results
    results_df = pd.concat(results)

    # Sort values but assign the result back to the variable
    results_df = results_df.sort_values(by="comparison")

    # Reset the index
    results_df = results_df.reset_index()

    # Rename index column to feature_id
    results_df = results_df.rename(columns={"index": "feature_id"})

    # Add taxonomic information
    daa_results_taxa_df = results_df.merge(
        taxonomy_table, how="left", left_on="feature_id", right_index=True
    )

    daa_results_taxa_fil_df = daa_results_taxa_df[daa_results_taxa_df["padj"] < 0.05]

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Save the results to a CSV file
    daa_results_taxa_df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

    # Save the filtered results to a CSV file
    daa_results_taxa_fil_df.to_csv(filtered_output_file, index=False)
    print(f"Results saved to {output_file} and {filtered_output_file}")

    return daa_results_taxa_df, daa_results_taxa_fil_df


def create_taxa_tables(
    filtered_asv_table_df: pd.DataFrame,
    filtered_taxonomy_df: pd.DataFrame,
    rank: str,
    output_file: str,
) -> pd.DataFrame:
    """
    Create taxa tables for species, genus, and family ranks and save them as CSV files.

    Args:
        filtered_asv_table_df (pd.DataFrame): Feature Filtered ASV table DataFrame.
        filtered_taxonomy_df (pd.DataFrame): Feature Filtered taxonomy DataFrame.
        rank (str): The taxonomic rank to process ('species', 'genus', or 'family').
        output_file (str, optional): Path to save the taxa table as CSV. If None, does not save.

    Returns:
        pd.DataFrame: Taxa table DataFrame.
    """

    # skip if output_file already exists
    if os.path.exists(output_file):
        print(f"✅ Taxa table for {rank} already exists at {output_file}, skipping.")

    # Initialize dictionary to store grouped data
    taxa_table_dict = {}

    # Group and sum ASV abundances by taxonomic rank
    for col in filtered_asv_table_df.columns:
        grouped_sums = (
            filtered_asv_table_df[col].groupby(filtered_taxonomy_df[rank]).sum()
        )
        taxa_table_dict[col] = grouped_sums

    # Convert to DataFrame
    taxa_table_df = pd.DataFrame(taxa_table_dict)

    # Save to CSV if output_file is provided
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    taxa_table_df.to_csv(output_file)
    print(f"Taxa table saved to: {output_file}")

    return taxa_table_df


def linda_daa_taxon(
    taxon_table_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    primary_variable: str,
    output_file: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the LINDA analysis for a single categorical comparison.

    Args:
        taxon_table (pd.DataFrame): The taxon count table.
        metadata (pd.DataFrame): The metadata.
        primary_variable (str): The column name for the categorical variable of interest.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: Full results and filtered significant results.
    """

    print(f"Running LINDA analysis for categorical variable: '{primary_variable}'")

    # Define path to filtered DAA results by replacing .csv extension with _filtered.csv
    filtered_output_file = output_file.replace(".csv", "_filtered.csv")

    # Skip if the filtered results file already exists
    if os.path.exists(filtered_output_file):
        print(
            f"Filtered results file {filtered_output_file} already exists. Skipping analysis."
        )
        return pd.read_csv(output_file), pd.read_csv(filtered_output_file)

    # Drop rows with missing values in the primary variable
    metadata_df = metadata_df.dropna(subset=[primary_variable])

    # Get unique groups in the primary variable
    unique_groups = metadata_df[primary_variable].unique()

    # Check if there are at least 2 unique groups
    if len(unique_groups) < 2:
        print(
            f"Error: Expected at least 2 unique groups in '{primary_variable}', found {len(unique_groups)}."
        )
        print("Skipping analysis due to insufficient groups.")
        return pd.DataFrame(), pd.DataFrame()

    # Generate all pairwise combinations of the unique groups
    pairwise_combinations = list(itertools.combinations(unique_groups, 2))

    # Create a list to store the results of the analysis
    results = []

    for group_a, group_b in pairwise_combinations:
        print(f"Processing pair: {group_a} vs {group_b}")
        # subset metadata to only include the two groups - CREATE EXPLICIT COPY
        metadata_subset_df = metadata_df[
            metadata_df[primary_variable].isin([group_a, group_b])
        ].copy()

        # Check if the subset has at least one sample for each group
        if metadata_subset_df[primary_variable].nunique() != 2:
            print(
                f"Skipping {group_a} vs {group_b}: insufficient data in metadata_subset"
            )
            continue

        # Subset OTU table base on metadata
        taxon_table_fil_df = taxon_table_df[metadata_subset_df.index]

        # Filter on prevalence
        taxon_table_fil_df = prevalance_filtering(taxon_table_fil_df, 10)

        # Skip if there are fewer than 4 samples (columns)
        if taxon_table_fil_df.shape[1] < 4:
            print(
                f"Error: only {taxon_table_fil_df.shape[1]} samples found (minimum 4 required)"
            )
            continue

        # Suppress FutureWarning from rpy2
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning)
            pandas2ri.activate()

            # Convert pandas dataframes to R dataframes
            r_taxon_table_df = pandas2ri.py2rpy(taxon_table_fil_df)
            r_metadata_subset_df = pandas2ri.py2rpy(metadata_subset_df)

        # Assigning converted dataframes to R's global environment
        r.assign("taxon_table_df", r_taxon_table_df)
        r.assign("metadata_df", r_metadata_subset_df)
        formula = f"~ {primary_variable}"
        r.assign("formula", formula)

        # Retrieve the factor levels of the primary variable in R
        levels_r = r(f"levels(as.factor(r_metadata${primary_variable}))")
        levels = list(levels_r)

        r_group_a, r_group_b = levels

        # Construct the path to the R script robustly
        current_script_dir = os.path.dirname(
            os.path.abspath(__file__)
        )  # Absolute path of the current script
        r_script_path = os.path.join(current_script_dir, "../R/daa_linda_taxa.R")

        r_script_path = os.path.abspath(r_script_path)  # Convert to an absolute path

        # Read the R code from the file
        with open(r_script_path, "r") as file:
            r_code = file.read()

        # Run the R code
        r(r_code)

        # Convert results to pandas dataframe - suppress warnings again
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning)
            linda_results = r["linda_results"]
            linda_results_df = pandas2ri.rpy2py(linda_results)

        # Add negative log10 p-value column
        linda_results_df["neg_log10_pvalue"] = -np.log10(linda_results_df["padj"])

        linda_results_df["direction"] = linda_results_df["log2FoldChange"].apply(
            lambda log2fc: f"Up in {r_group_b}" if log2fc > 0 else f"Up in {r_group_a}"
        )

        linda_results_df["comparison"] = f"{r_group_a} vs {r_group_b}"
        results.append(linda_results_df)

    # Combine all results
    results_df = pd.concat(results)

    # Sort values
    results_df = results_df.sort_values(by="padj")

    # Reset the index
    results_df = results_df.reset_index()

    # Rename index column to feature_id
    results_df = results_df.rename(columns={"index": "feature_id"})

    # Filter for significant results
    results_filtered_df = results_df[results_df["padj"] < 0.05]

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    # Save the results to a CSV file
    results_df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")
    # Save the filtered results to a CSV file
    results_filtered_df.to_csv(filtered_output_file, index=False)
    print(f"Filtered results saved to {filtered_output_file}")

    return results_df, results_filtered_df


def linda_daa_picrust2(
    function_table_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    primary_variable: str,
    output_file: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run the LINDA analysis for a single categorical comparison.

    Args:
        taxon_table (pd.DataFrame): The taxon count table.
        metadata (pd.DataFrame): The metadata.
        primary_variable (str): The column name for the categorical variable of interest.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: Full results and filtered significant results.
    """

    print(f"Running LINDA analysis for categorical variable: '{primary_variable}'")

    # Define path to filtered DAA results by replacing .csv extension with _filtered.csv
    filtered_output_file = output_file.replace(".csv", "_filtered.csv")

    # Skip if the filtered results file already exists
    if os.path.exists(filtered_output_file):
        print(
            f"Filtered results file {filtered_output_file} already exists. Skipping analysis."
        )
        return pd.read_csv(output_file), pd.read_csv(filtered_output_file)

    # Drop rows with missing values in the primary variable
    metadata_df = metadata_df.dropna(subset=[primary_variable])

    # Get unique groups in the primary variable
    unique_groups = metadata_df[primary_variable].unique()

    # Check if there are at least 2 unique groups
    if len(unique_groups) < 2:
        print(
            f"Error: Expected at least 2 unique groups in '{primary_variable}', found {len(unique_groups)}."
        )
        print("Skipping analysis due to insufficient groups.")
        return pd.DataFrame(), pd.DataFrame()

    # Generate all pairwise combinations of the unique groups
    pairwise_combinations = list(itertools.combinations(unique_groups, 2))

    # Create a list to store the results of the analysis
    results = []

    for group_a, group_b in pairwise_combinations:
        print(f"Processing pair: {group_a} vs {group_b}")
        # subset metadata to only include the two groups - CREATE EXPLICIT COPY
        metadata_subset_df = metadata_df[
            metadata_df[primary_variable].isin([group_a, group_b])
        ].copy()

        # Check if the subset has at least one sample for each group
        if metadata_subset_df[primary_variable].nunique() != 2:
            print(
                f"Skipping {group_a} vs {group_b}: insufficient data in metadata_subset"
            )
            continue

        # Subset OTU table base on metadata
        function_table_fil_df = function_table_df[metadata_subset_df.index]

        # Filter on prevalence
        function_table_fil_df = prevalance_filtering(function_table_fil_df, 10)

        # Skip if there are fewer than 4 samples (columns)
        if function_table_fil_df.shape[1] < 4:
            print(
                f"Error: only {function_table_fil_df.shape[1]} samples found (minimum 4 required)"
            )
            continue

        # Suppress FutureWarning from rpy2
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning)
            pandas2ri.activate()

            # Convert pandas dataframes to R dataframes
            r_function_table_df = pandas2ri.py2rpy(function_table_fil_df)
            r_metadata_subset_df = pandas2ri.py2rpy(metadata_subset_df)

        # Assigning converted dataframes to R's global environment
        r.assign("taxon_table_df", r_function_table_df)
        r.assign("metadata_df", r_metadata_subset_df)
        formula = f"~ {primary_variable}"
        r.assign("formula", formula)

        # Retrieve the factor levels of the primary variable in R
        levels_r = r(f"levels(as.factor(metadata_df${primary_variable}))")
        levels = list(levels_r)

        r_group_a, r_group_b = levels

        # Construct the path to the R script robustly
        current_script_dir = os.path.dirname(
            os.path.abspath(__file__)
        )  # Absolute path of the current script
        r_script_path = os.path.join(current_script_dir, "../R/daa_linda_taxa.R")

        r_script_path = os.path.abspath(r_script_path)  # Convert to an absolute path

        # Read the R code from the file
        with open(r_script_path, "r") as file:
            r_code = file.read()

        # Run the R code
        r(r_code)

        # Convert results to pandas dataframe - suppress warnings again
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning)
            linda_results = r["linda_results"]
            linda_results_df = pandas2ri.rpy2py(linda_results)

        # Add negative log10 p-value column
        linda_results_df["neg_log10_pvalue"] = -np.log10(linda_results_df["padj"])

        linda_results_df["direction"] = linda_results_df["log2FoldChange"].apply(
            lambda log2fc: f"Up in {r_group_b}" if log2fc > 0 else f"Up in {r_group_a}"
        )

        linda_results_df["comparison"] = f"{r_group_a} vs {r_group_b}"
        results.append(linda_results_df)

    # Combine all results
    results_df = pd.concat(results)

    # Sort values
    results_df = results_df.sort_values(by="padj")

    # Reset the index
    results_df = results_df.reset_index()

    # Rename index column to feature_id
    results_df = results_df.rename(columns={"index": "feature_id"})

    # Filter for significant results
    results_filtered_df = results_df[results_df["padj"] < 0.05]

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    # Save the results to a CSV file
    results_df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")
    # Save the filtered results to a CSV file
    results_filtered_df.to_csv(filtered_output_file, index=False)
    print(f"Filtered results saved to {filtered_output_file}")

    return results_df, results_filtered_df
