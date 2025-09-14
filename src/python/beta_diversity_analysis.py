import os

import biom  # type: ignore
import pandas as pd
import unifrac  # type: ignore
from biom.table import Table  # type: ignore
from gemelli.rpca import phylogenetic_rpca, rpca  # type: ignore
from rpy2.robjects import pandas2ri, r  # type: ignore
from skbio import DistanceMatrix, OrdinationResults, TreeNode  # type: ignore
from skbio.diversity import beta_diversity  # type: ignore

from python.processing import reformat_taxonomy


def calculate_jaccard_distance(
    rarefied_table_df: pd.DataFrame, output_file: str, qimme2_format: bool = True
) -> tuple[pd.DataFrame, DistanceMatrix]:
    """
    Perform Jaccard distance calculation.

    Parameters:
    otu_table (pd.DataFrame): OTU table with samples as rows and OTUs as columns.
    sample_ids (list): List of sample IDs.

    Returns:
    pd.DataFrame: Jaccard distance matrix.
    """

    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    if qimme2_format:
        # QIIME2 formats the OTU table with samples as rows and OTUs as columns
        rarefied_table_df = rarefied_table_df.T

    # Calculate Jaccard distance matrix
    jaccard_dm = beta_diversity(
        metric="jaccard",
        counts=rarefied_table_df.values,
        ids=rarefied_table_df.index.tolist(),
    )

    # Convert DistanceMatrix to DataFrame
    jaccard_df = jaccard_dm.to_data_frame()

    # Write the Jaccard distance matrix
    jaccard_df.to_csv(output_file)
    print(f"Jaccard distance matrix saved to {output_file}")

    return jaccard_df, jaccard_dm


def calculate_bray_curtis_distance(
    rarefied_table_df: pd.DataFrame,
    output_file: str,
) -> tuple[pd.DataFrame, DistanceMatrix]:
    """Perform Bray-Curtis distance calculation on the rarefied OTU table.

    Args:
        rarefied_otu_table (pd.DataFrame): The rarefied OTU table where rows are taxa and columns are samples.

    Returns:
        pd.DataFrame: A DataFrame representing the Bray-Curtis distance matrix.

    """

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Transpose the OTU table to have samples as rows and taxa as columns (as required by Bray-Curtis)
    otu_data = rarefied_table_df.T.to_numpy()

    # Perform Bray-Curtis calculation
    bray_curtis_dm = beta_diversity(
        metric="braycurtis", counts=otu_data, ids=rarefied_table_df.columns
    )

    # Convert to data-frame
    bray_curtis_df = bray_curtis_dm.to_data_frame()

    # Write the Bray-Curtis distance matrix
    bray_curtis_df.to_csv(output_file)
    print(f"Bray-Curtis distance matrix saved to {output_file}")

    return bray_curtis_df, bray_curtis_dm


def perform_generalized_unifrac_distance(
    biom_table_path: str,
    rooted_tree_newick_path: str,
    output_file: str,
    alpha: float = 0.5,
) -> tuple[pd.DataFrame, DistanceMatrix]:
    """Calculate Unifrac distance matrix from a BIOM table and a rooted tree.

    Args:
        biom_table_path (str): Path to the BIOM file (must be in BIOM 2.1 format).
        rooted_tree_newick_path (str): Path to the rooted tree in Newick format.
        output_dir (str): Directory to save the output distance matrix.
        alpha (float, optional): Generalized UniFrac alpha parameter (default 0.5).

    Returns:
        tuple[pd.DataFrame, DistanceMatrix]: A DataFrame and DistanceMatrix representing the UniFrac distance matrix.
    """

    # Compute Generalized UniFrac distance matrix
    unifrac_distance_matrix = unifrac.generalized(
        table=biom_table_path,
        phylogeny=rooted_tree_newick_path,
        alpha=alpha,
    )

    # Convert to DataFrame
    unifrac_distance_df = unifrac_distance_matrix.to_data_frame()

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Write the UniFrac distance matrix to CSV
    unifrac_distance_df.to_csv(output_file, index=False)
    print(f"Generalized UniFrac distance matrix saved to {output_file}")

    return unifrac_distance_df, unifrac_distance_matrix


def perform_permanova_with_vegan(
    metadata: pd.DataFrame,
    distance_df: pd.DataFrame,
    primary_variable: str,
    strata: str = "False",
) -> tuple[dict, pd.DataFrame]:
    """Performs PERMANOVA analysis and tests homogeneity of dispersion using betadisper.

    Args:
        metadata (pd.DataFrame): A metadata table containing sample information.
        distance_matrix (pd.DataFrame): A distance matrix representing the distances between samples.
        primary_variable (str): The primary variable/column in the metadata to test.

    Returns:
        tuple[dict, pd.DataFrame]: A dictionary containing PERMANOVA and betadisper results and the results DataFrame.
    """

    # Activate pandas-to-R conversion
    pandas2ri.activate()

    # Convert pandas DataFrame to R dataframe
    r_metadata = pandas2ri.py2rpy(metadata)
    r_distance_matrix = pandas2ri.py2rpy(distance_df)

    # Assign the metadata and distance matrix to the R environment
    r.assign("distance_matrix", r_distance_matrix)
    r.assign("metadata", r_metadata)
    r.assign("primary_variable", primary_variable)
    r.assign("strata", strata)

    # Construct the path to the R script robustly
    current_script_dir = os.path.dirname(
        os.path.abspath(__file__)
    )  # Absolute path of the current script
    r_script_path = os.path.join(current_script_dir, "../R/permanova_vegan.R")

    r_script_path = os.path.abspath(r_script_path)  # Convert to an absolute path

    # Read the R code from the file
    with open(r_script_path, "r") as file:
        r_code = file.read()

    # Run the R code
    r(r_code)

    # Retrieve the results data frame from R
    results_df = r["results_df"]
    results_df = pandas2ri.rpy2py(results_df)

    # Convert the R data frame to a Python dictionary for easy use
    results_dict = results_df.set_index("metric")["value"].to_dict()

    print(f"PERMANOVA and betadisper results: {results_dict}")
    return results_dict


def perform_rpca(
    asv_table_df: pd.DataFrame,
    output_file: str,
) -> tuple[OrdinationResults, DistanceMatrix]:
    """Perform Robust Aitchison PCA (RPCA) on BIOM table.

    Args:
        table (Table): BIOM table

    Returns:
        tuple[OrdinationResults, DistanceMatrix]: coordinates and distance matrix
    """

    # Convert dataframe to biom Table
    biom_table = biom.Table(
        asv_table_df.values,
        observation_ids=asv_table_df.index.astype(str),
        sample_ids=asv_table_df.columns.astype(str),
    )

    # Preform RPCA
    ordination, robust_aitchison_dm = rpca(biom_table, min_sample_count=500)

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Write outputs
    ordination.write(
        output_file.replace("distance_matrix.csv", "ordination.txt"),
        format="ordination",
    )

    # Convert distance matrix to DataFrame
    distance_df = robust_aitchison_dm.to_data_frame()

    # Write out distance DataFrame to CSV
    distance_df.to_csv(output_file)

    return ordination, distance_df


def perform_phylogenetic_rpca(
    asv_table_df: pd.DataFrame,
    rooted_tree_newick_path: str,
    taxonomy_df: pd.DataFrame,
    output_file: str,
) -> tuple[OrdinationResults, pd.DataFrame, TreeNode, pd.DataFrame, pd.DataFrame]:
    """Perform Phylogenetic RPCA on ASV table, tree, and taxonomy.

    Args:
        asv_table_df (pd.DataFrame): Feature table (features x samples).
        rooted_tree_newick_path (str): Path to the rooted tree in Newick format.
        taxonomy_df (pd.DataFrame): Taxonomy DataFrame (required).
        output_file (str): Path to write out pandas dataframe as csv.

    Returns:
        tuple: ordination results, distance DataFrame, pruned tree, filtered table, filtered taxonomy.
    """
    # Convert DataFrame to BIOM Table
    biom_table = Table(
        asv_table_df.values,
        observation_ids=asv_table_df.index.astype(str),
        sample_ids=asv_table_df.columns.astype(str),
    )

    taxonomy_df = taxonomy_df[["Taxon"]]  # Ensure only the 'Taxon' column is present

    # Run Phylogenetic RPCA with taxonomy

    (
        ordination_results,
        phylo_rclr_distance_matrix,
        pruned_tree,
        filtered_table,
        filtered_taxonomy,
    ) = phylogenetic_rpca(
        table=biom_table,
        phylogeny=rooted_tree_newick_path,
        taxonomy=taxonomy_df,
        min_sample_count=500,
    )

    filtered_taxonomy = reformat_taxonomy(taxonomy_df=filtered_taxonomy)

    # Convert distance matrix to DataFrame
    distance_df = phylo_rclr_distance_matrix.to_data_frame()

    # Write outputs
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Write distance DataFrame to CSV
    distance_df.to_csv(output_file, index=False)

    ordination_results.write(
        output_file.replace("distance_matrix.csv", "ordination.txt"),
        format="ordination",
    )

    return (
        ordination_results,
        distance_df,
        pruned_tree,
        filtered_table,
        filtered_taxonomy,
    )


def calculate_unweighted_unifrac_distance(
    biom_table_path: str,
    rooted_tree_newick_path: str,
    output_file: str,
) -> tuple[pd.DataFrame, DistanceMatrix]:
    """Calculate Unweighted Unifrac distance matrix from a BIOM table and a rooted tree.

    Args:
        biom_table_path (str): Path to the BIOM file (must be in BIOM 2.1 format).
        rooted_tree_newick_path (str): Path to the rooted tree in Newick format.
        output_dir (str): Directory to save the output distance matrix.

    Returns:
        tuple[pd.DataFrame, DistanceMatrix]: A DataFrame and DistanceMatrix representing the UniFrac distance matrix.
    """

    # Compute Unweighted UniFrac distance matrix
    unifrac_distance_matrix = unifrac.unweighted(
        table=biom_table_path,
        phylogeny=rooted_tree_newick_path,
    )

    # Convert to DataFrame
    unifrac_distance_df = unifrac_distance_matrix.to_data_frame()

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Write the UniFrac distance matrix to CSV
    unifrac_distance_df.to_csv(output_file, index=False)
    print(f"Unweighted UniFrac distance matrix saved to {output_file}")

    return unifrac_distance_df, unifrac_distance_matrix


def perform_pairwise_comparisons_permanova(
    metadata: pd.DataFrame,
    distance_df: pd.DataFrame,
    primary_variable: str,
    output_file: str,
    all: bool = False,
) -> None:
    """Performs PERMANOVA analysis and tests homogeneity of dispersion using betadisper.

    Args:
        metadata (pd.DataFrame): A metadata table containing sample information.
        distance_matrix (pd.DataFrame): A distance matrix representing the distances between samples.

    Returns:
        Dict[str, Any]: A dictionary containing PERMANOVA and betadisper results, including pseudo-F statistic, p-value, R²,
                        dispersion F-statistic, and dispersion p-value.
    """

    if all:
        strata = "True"
    else:
        strata = "False"

    # Activate pandas-to-R conversion
    pandas2ri.activate()

    # Convert pandas DataFrame to R dataframe
    r_metadata = pandas2ri.py2rpy(metadata)
    r_distance_matrix = pandas2ri.py2rpy(distance_df)

    # Assign the metadata and distance matrix to the R environment
    r.assign("distance_matrix", r_distance_matrix)
    r.assign("metadata", r_metadata)
    r.assign("primary_variable", primary_variable)
    r.assign("strata", strata)
    r.assign("output_file", output_file)

    # Construct the path
    # Get current path
    current_script_dir = os.path.dirname(os.path.abspath(__file__))
    # Construct the path to the R script
    r_script_path = os.path.join(
        current_script_dir, "../R/pairwise_comparisons_adonis.R"
    )
    # Ensure the path is absolute
    r_script_path = os.path.abspath(r_script_path)

    # Read the R code from the provided script
    with open(r_script_path, "r") as file:
        r_code = file.read()

    # Execute the R code
    r(r_code)
