import os
import warnings
from itertools import combinations
from typing import Optional

import biom  # type: ignore
import numpy as np
import pandas as pd
import scipy.stats  # type: ignore
import unifrac  # type: ignore
from biom.table import Table  # type: ignore
from gemelli.ctf import ctf  # type: ignore
from gemelli.rpca import phylogenetic_rpca, rpca  # type: ignore
from qiime2 import Artifact  # type: ignore
from qiime2.plugins.diversity.pipelines import alpha_phylogenetic  # type: ignore
from qiime2.plugins.phylogeny import pipelines as phylo_pipelines  # type: ignore
from rpy2 import robjects as r  # type: ignore
from rpy2.robjects import pandas2ri  # type: ignore
from skbio import DistanceMatrix, OrdinationResults, TreeNode  # type: ignore
from skbio.diversity import alpha, beta_diversity  # type: ignore
from sklearn.manifold import MDS  # type: ignore

from python.processing import reformat_taxonomy


def depth_filtering(table: pd.DataFrame, depth_threshold: float) -> pd.DataFrame:
    """Remove samples from feature table with read depth below threshold."""
    sample_depths = table.sum(axis=0)

    # Remove samples with read depth below threshold
    return table.loc[:, sample_depths >= depth_threshold]


def generate_phylogenetic_tree_with_qiime2(input_fasta: str, threads, output_dir: str):
    """
    Generates a phylogenetic tree from an input FASTA file using QIIME2 MAFFT and FastTree.

    Returns:
    - rooted_tree (Artifact): The rooted tree QIIME2 artifact object.
    - unrooted_tree (Artifact): The unrooted tree QIIME2 artifact object.
    - rooted_tree_newick (str): Path to the exported rooted tree in Newick format.
    """
    rooted_tree_file = os.path.join(output_dir, "rooted_tree.qza")
    unrooted_tree_file = os.path.join(output_dir, "unrooted_tree.qza")
    rooted_tree_newick_path = os.path.join(output_dir, "tree.nwk")

    # Skip if output directory already exists and read in trees
    if os.path.exists(rooted_tree_file) and os.path.exists(rooted_tree_newick_path):
        print(f"Output files already exist in {output_dir}. Skipping tree generation.")
        print("Reading in existing trees...")
        rooted_tree = Artifact.load(rooted_tree_file)
        unrooted_tree = Artifact.load(unrooted_tree_file)
        return rooted_tree, unrooted_tree, rooted_tree_newick_path

    print("Starting to generate phylogenetic tree...")
    os.makedirs(output_dir, exist_ok=True)

    rep_seq_artifact = Artifact.import_data(
        type="FeatureData[Sequence]", view=input_fasta
    )

    print("✅ Finished reading in sequences.")

    action_results = phylo_pipelines.align_to_tree_mafft_fasttree(
        sequences=rep_seq_artifact, n_threads=threads
    )

    rooted_tree = action_results.rooted_tree
    unrooted_tree = action_results.tree

    rooted_tree.save(rooted_tree_file)
    unrooted_tree.save(unrooted_tree_file)
    print(f"Rooted tree saved to {rooted_tree_file}")
    print(f"Unrooted tree saved to {unrooted_tree_file}")

    # Export rooted tree to Newick format
    rooted_tree.export_data(output_dir)
    print(f"Rooted tree exported to Newick format at {output_dir}")

    return rooted_tree, unrooted_tree, rooted_tree_newick_path


def rarefy_df(
    df: pd.DataFrame,
    sampling_depth: Optional[int] = None,
    with_replacement: bool = False,
    random_seed: int = 1088,
) -> pd.DataFrame:
    """
    Rarefy a pandas DataFrame to a specified sampling depth using the BIOM format.
    This function converts a DataFrame to a BIOM table, performs rarefaction,
    and returns the rarefied BIOM table along with the DataFrame representation.
    The DataFrame should have samples as columns and features as rows.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with samples as columns and features as rows.
        The DataFrame should contain non-negative integer counts.

    sampling_depth : int, optional
        Rarefaction depth. If None, automatically set to the minimum sample sum.

    with_replacement : bool, default False
        Perform rarefaction sampling with or without replacement.

    random_seed : int, default 1088
        Seed for random number generation to ensure reproducibility.

    Returns
    -------
    pd.DataFrame
        The rarefied DataFrame.
    """

    if sampling_depth is None:
        # Get the minimum sample sum for rarefaction depth from the DataFrame
        sampling_depth = int(min(df.sum(axis=0)))
    else:
        sampling_depth = int(sampling_depth)

    print(f"Sampling depth set to {sampling_depth}.")

    # Convert dataframe to biom Table
    biom_table = biom.Table(
        df.values,
        observation_ids=df.index.astype(str),
        sample_ids=df.columns.astype(str),
    )

    biom_table = biom_table.filter(
        lambda v, i, m: v.sum() >= sampling_depth, inplace=False, axis="sample"
    )
    rarefied_biom_table = biom_table.subsample(
        sampling_depth,
        axis="sample",
        by_id=False,
        with_replacement=with_replacement,
        seed=random_seed,
    )

    # Convert back to DataFrame for consistency
    rarefied_df = pd.DataFrame(
        rarefied_biom_table.matrix_data.toarray(),
        index=rarefied_biom_table.ids(axis="observation"),
        columns=rarefied_biom_table.ids(axis="sample"),
    )

    print(
        f"Rarefied table shape: {rarefied_df.shape}, "
        f"with {rarefied_df.shape[0]} features and {rarefied_df.shape[1]} samples."
    )
    print(f"Sampling depth used for rarefaction: {sampling_depth}.")

    return rarefied_df


def convert_dataframe_to_qiime2_artifact(df: pd.DataFrame) -> Artifact:
    """
    Convert a pandas DataFrame to a QIIME2 Artifact of type 'FeatureTable[Frequency]'.
    """
    biom_table = biom.Table(
        df.values,
        observation_ids=df.index.astype(str),
        sample_ids=df.columns.astype(str),
    )
    artifact = Artifact.import_data("FeatureTable[Frequency]", biom_table)
    return artifact


def calculate_phylo_alpha_diversity_qiime2(
    normalized_table: Artifact, rooted_tree: Artifact, metric: str = "faith_pd"
) -> pd.DataFrame:
    """
    Calculate alpha diversity using a phylogenetic metric and return a DataFrame.

    Parameters:
    normalized_table (pd.DataFrame): The rarefied/normalized feature table.
    rooted_tree: The phylogenetic tree associated with the dataset.
    metric (str): The diversity metric to use (default is 'faith_pd').

    Returns:
    pd.DataFrame: A DataFrame containing alpha diversity values.
    """
    ap_results = alpha_phylogenetic(
        table=normalized_table, phylogeny=rooted_tree, metric=metric
    )
    ap_numbers = ap_results.alpha_diversity
    alpha_diversity_series = ap_numbers.view(pd.Series)
    alpha_diversity_df = alpha_diversity_series.to_frame(
        name="Faith's phylogenetic diversity"
    )

    # Reset index and rename columns
    alpha_diversity_df.reset_index(inplace=True)
    alpha_diversity_df.rename(columns={"index": "sample_id"}, inplace=True)

    return alpha_diversity_df


def calculate_alpha_diversity(
    normalized_table_df: pd.DataFrame,
    normalized_table_qza: Artifact,
    rooted_tree: Artifact,
    output_file: str,
) -> pd.DataFrame:
    """Calculates alpha diversity metrics for each sample in the normalized table.

    Args:
        normalized_table (pd.DataFrame): The normalized OTU table where rows are taxa and columns are samples.

    Returns:
        pd.DataFrame: A DataFrame containing alpha diversity metrics for each sample.
    """

    # Calculate alpha diversity metrics
    print("Calculating alpha diversity metrics...")
    alpha_diversity_metrics: dict[str, list] = {
        "sample_id": [],
        "Observed Features": [],
        "Simpson's dominance index": [],
        "Simpson's diversity index": [],
    }

    for sample in normalized_table_df.columns:
        # Calculate observed OTUs
        observed_otus = alpha.observed_otus(normalized_table_df[sample])

        # Calculate Simpson's dominance index
        dominance = alpha.dominance(normalized_table_df[sample])

        # Calculate Simpson diversity index
        simpson = alpha.simpson(normalized_table_df[sample])

        alpha_diversity_metrics["sample_id"].append(sample)
        alpha_diversity_metrics["Observed Features"].append(observed_otus)
        alpha_diversity_metrics["Simpson's dominance index"].append(dominance)
        alpha_diversity_metrics["Simpson's diversity index"].append(simpson)

    alpha_diversity_df = pd.DataFrame(alpha_diversity_metrics)

    # Calculate Faith's phylogenetic diversity if a tree is provided
    phylo_alpha_metrics_df = calculate_phylo_alpha_diversity_qiime2(
        normalized_table_qza, rooted_tree, metric="faith_pd"
    )

    alpha_diversity_merged_df = phylo_alpha_metrics_df.merge(
        alpha_diversity_df, left_on="sample_id", right_on="sample_id"
    )

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    # Save the alpha diversity metrics to a CSV file
    alpha_diversity_merged_df.to_csv(output_file, index=False)
    print(f"Alpha diversity metrics saved to {output_file}")

    return alpha_diversity_merged_df


def pairwise_alpha_diversity_calculations(
    alpha_diversity_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    primary_variable: str,
    output_file: str,
) -> pd.DataFrame:
    """Perform pairwise comparisons of alpha diversity between groups.

    Args:
        alpha_diversity_df (pd.DataFrame): DataFrame containing alpha diversity metrics.
        metadata_df (pd.DataFrame): DataFrame containing sample metadata.
        primary_variable (str): The column in metadata_df to group by.
        output_file (str): Path to write out pandas dataframe as csv.

    Returns:
        pd.DataFrame: DataFrame containing pairwise comparison results.
    """

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Initialize results list to store comparison results
    results_list = []

    # Get groups from metadata
    groups = metadata_df[primary_variable].unique()

    for g1, g2 in combinations(groups, 2):
        print(g1, g2)
        comparison = f"{g1}_vs_{g2}"

        # Subset metadata by g1
        meta_g1 = metadata_df[metadata_df[primary_variable] == g1]
        meta_g2 = metadata_df[metadata_df[primary_variable] == g2]

        # Find the intersection of patients_id in both groups
        common_patients = set(meta_g1["patient_id"]).intersection(
            set(meta_g2["patient_id"])
        )

        if not common_patients:
            print(
                f"No common patients found between groups {g1} and {g2}. Skipping comparison."
            )
            continue

        # Subset metadata to only include common patients
        meta_data_subset_df = metadata_df[
            metadata_df["patient_id"].isin(common_patients)
        ]

        # Merge alpha diversity with subset metadata
        merged_df = alpha_diversity_df.merge(
            meta_data_subset_df, left_on="sample_id", right_on="sample_id"
        )

        # Perform pairwise wilcoxon test for each alpha diversity metric

        for metric in [
            "Faith's phylogenetic diversity",
            "Observed Features",
            "Simpson's dominance index",
            "Simpson's diversity index",
        ]:
            group1_values = merged_df[merged_df[primary_variable] == g1][metric]
            group2_values = merged_df[merged_df[primary_variable] == g2][metric]

            # Perform Wilcoxon rank-sum test (Mann-Whitney U test)
            stat, p_value = scipy.stats.mannwhitneyu(
                group1_values, group2_values, alternative="two-sided"
            )

            # Diffenece in medians
            median_diff = group1_values.median() - group2_values.median()

            if median_diff > 0:
                relative_change = "up in " + str(g1)
            elif median_diff < 0:
                relative_change = "up in " + str(g2)
            else:
                relative_change = "no change"

            # append results to a list
            results_list.append(
                {
                    "group1": g1,
                    "group2": g2,
                    "comparison": comparison,
                    "metric": metric,
                    "relative_change": relative_change,
                    "U-statistic": stat,
                    "p-value": p_value,
                    "median_diff": median_diff,
                }
            )

    # Convert results list to DataFrame
    results_df = pd.DataFrame(results_list)

    # Add a column for significance based on p-value via stars
    results_df["significance"] = results_df["p-value"].apply(
        lambda p: "****"
        if p < 0.0001
        else "***"
        if p < 0.001
        else "**"
        if p < 0.01
        else "*"
        if p < 0.05
        else "ns"
    )

    # Write results DataFrame to CSV
    results_df.to_csv(output_file, index=False)

    return results_df


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


def perform_phylo_rpca(
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

    taxonomy_df = taxonomy_df[["Taxon", "Confidence"]]

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


def perform_bray_curtis(
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
    bray_curtis_results = beta_diversity(
        metric="braycurtis", counts=otu_data, ids=rarefied_table_df.columns
    )

    # Convert to a DistanceMatrix object for further downstream analysis
    bray_curtis_dm = DistanceMatrix(
        bray_curtis_results.data, ids=bray_curtis_results.ids
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
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
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
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        results_df = r["results_df"]
        results_df = pandas2ri.rpy2py(results_df)

    # Convert the R data frame to a Python dictionary for easy use
    results_dict = results_df.set_index("metric")["value"].to_dict()

    print(f"PERMANOVA and betadisper results: {results_dict}")
    return results_dict


def perform_nmds(distance_matrix: DistanceMatrix) -> tuple[pd.DataFrame, float]:
    """Perform NMDS (Non-metric Multidimensional Scaling) on a Bray-Curtis distance matrix.

    Args:
        distance_matrix (DistanceMatrix): A distance matrix representing the Bray-Curtis distances between samples.

    Returns:
        tuple[pd.DataFrame, float]: A DataFrame containing sample coordinates and the stress value.
    """
    # Convert skbio DistanceMatrix to a numpy array for use with sklearn
    dist_array = distance_matrix.data

    # Use MDS from sklearn for NMDS with precomputed distances
    mds = MDS(
        n_components=2,
        metric=False,
        dissimilarity="precomputed",
        random_state=1835,
        n_init=10,
    )
    coords = mds.fit_transform(dist_array)

    # Return coordinates as DataFrame and the stress value
    coordinates_df = pd.DataFrame(
        coords, index=distance_matrix.ids, columns=["NMDS1", "NMDS2"]
    )

    normalized_stress = np.sqrt(mds.stress_ / np.sum(dist_array**2))

    return coordinates_df, normalized_stress


def perform_ctf(
    asv_count_table_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    individual_id_column: str,
    state_column: str,
    output_dir: str,
) -> tuple[OrdinationResults, pd.DataFrame]:
    """Perform Compositional Tensor Factorization (CTF) on BIOM table.

    Args:
        asv_count_table_biom: BIOM table (feature table)
        metadata_df (pd.DataFrame): sample metadata
        individual_id_column (str): subject identifier column
        state_column (str): state/timepoint column
        output_dir (str): not currently used

    Returns:
        tuple[
            DistanceMatrix,   # distance matrix object (skbio)
            pd.DataFrame,     # distance matrix as DataFrame
            pd.DataFrame,     # sample coordinates
            pd.DataFrame      # feature coordinates
        ]
    """

    # Filter asv_count_table_df for depth
    asv_count_table_df = depth_filtering(asv_count_table_df, depth_threshold=1000)

    # Convert dataframe to biom Table
    biom_table = biom.Table(
        asv_count_table_df.values,
        observation_ids=asv_count_table_df.index.astype(str),
        sample_ids=asv_count_table_df.columns.astype(str),
    )

    # Run CTF
    ordination_results, feature_ordination, distance_matrix, sample_df, feature_df = (
        ctf(
            biom_table.copy(),
            sample_metadata=metadata_df,
            individual_id_column=individual_id_column,
            state_column=state_column,
            n_components=2,
            min_sample_count=500,
        )
    )
    print("CTF complete.")

    # Convert distance matrix to DataFrame
    distance_df = distance_matrix.to_data_frame()

    # Write output to output_dir
    # make output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    distance_matrix.write(os.path.join(output_dir, "ctf_distance_matrix.txt"))
    ordination_results.write(
        os.path.join(output_dir, "ctf_sample_ordination.txt"), format="ordination"
    )
    feature_ordination.write(
        os.path.join(output_dir, "ctf_feature_ordination.txt"), format="ordination"
    )
    sample_df.to_csv(os.path.join(output_dir, "ctf_sample_coordinates.csv"))
    feature_df.to_csv(os.path.join(output_dir, "ctf_feature_coordinates.csv"))

    return ordination_results, distance_df


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
