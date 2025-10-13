import os
from itertools import combinations
from typing import Optional

import biom  # type: ignore
import numpy as np
import pandas as pd
import scipy.stats  # type: ignore

# Disable CUDA to avoid GPU initialization errors
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["NUMBA_DISABLE_CUDA"] = "1"

from gemelli.ctf import ctf  # type: ignore
from qiime2 import Artifact  # type: ignore
from qiime2.plugins.diversity.pipelines import alpha_phylogenetic  # type: ignore
from qiime2.plugins.phylogeny import pipelines as phylo_pipelines  # type: ignore
from skbio import DistanceMatrix, OrdinationResults  # type: ignore
from skbio.diversity import alpha  # type: ignore
from sklearn.manifold import MDS  # type: ignore


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


def paired_alpha_diversity_calculations(
    alpha_diversity_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    primary_variable: str,
    output_file: str,
) -> pd.DataFrame:
    """
    Perform pairwise comparisons of alpha diversity between groups using paired Wilcoxon signed-rank tests.

    Args:
        alpha_diversity_df (pd.DataFrame): DataFrame containing alpha diversity metrics (indexed by sample_id).
        metadata_df (pd.DataFrame): DataFrame containing sample metadata (including 'sample_id' and 'patient_id').
        primary_variable (str): The column in metadata_df to group by (e.g., condition, timepoint).
        output_file (str): Path to write the output CSV file.

    Returns:
        pd.DataFrame: DataFrame containing pairwise comparison results with z-score scaled within-patient differences.
    """

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    results_list = []

    # Get unique groups for comparisons
    groups = metadata_df[primary_variable].dropna().unique()

    # Loop over all pairwise group comparisons
    for g1, g2 in combinations(groups, 2):
        print(f"Comparing {g1} vs {g2}")

        # Get patients that appear in both groups
        meta_g1 = metadata_df[metadata_df[primary_variable] == g1]
        meta_g2 = metadata_df[metadata_df[primary_variable] == g2]
        common_patients = set(meta_g1["patient_id"]).intersection(meta_g2["patient_id"])

        if not common_patients:
            print(f"No common patients between {g1} and {g2}. Skipping.")
            continue

        # Subset to common patients only
        subset_meta = metadata_df[metadata_df["patient_id"].isin(common_patients)]

        # Merge alpha diversity data with metadata
        merged_df = alpha_diversity_df.merge(subset_meta, on="sample_id")

        # Define metrics to test
        metrics = [
            "Faith's phylogenetic diversity",
            "Observed Features",
            "Simpson's dominance index",
            "Simpson's diversity index",
        ]

        for metric in metrics:
            # Pivot so each row = patient, columns = group values
            paired_df = merged_df.pivot_table(
                index="patient_id", columns=primary_variable, values=metric
            ).dropna(subset=[g1, g2])

            if paired_df.shape[0] < 1:
                print(f"No paired data for {metric} in {g1} vs {g2}. Skipping.")
                continue

            # Extract aligned values
            x = paired_df[g1]
            y = paired_df[g2]

            # Paired Wilcoxon signed-rank test
            stat, p_value = scipy.stats.wilcoxon(x, y, alternative="two-sided")

            # Compute within-patient differences
            differences = x - y
            raw_median_diff = differences.median()
            pooled_std = differences.std()
            zscore_diff = raw_median_diff / pooled_std if pooled_std != 0 else 0

            # Direction of change
            if raw_median_diff > 0:
                relative_change = f"up in {g1}"
            elif raw_median_diff < 0:
                relative_change = f"up in {g2}"
            else:
                relative_change = "no change"

            # Store results
            results_list.append(
                {
                    "group1": g1,
                    "group2": g2,
                    "comparison": f"{g1}_vs_{g2}",
                    "metric": metric,
                    "relative_change": relative_change,
                    "Wilcoxon_statistic": stat,
                    "p-value": p_value,
                    "median_diff": raw_median_diff,
                    "zscore_diff": zscore_diff,
                }
            )

    # Convert results list to DataFrame
    results_df = pd.DataFrame(results_list)

    # Significance stars
    def get_significance(p):
        if p < 0.0001:
            return "****"
        elif p < 0.001:
            return "***"
        elif p < 0.01:
            return "**"
        elif p < 0.05:
            return "*"
        else:
            return "ns"

    results_df["significance"] = results_df["p-value"].apply(get_significance)

    # Write to CSV
    results_df.to_csv(output_file, index=False)
    print(f"\nResults saved to {output_file}")

    return results_df


def unpaired_alpha_diversity_calculations(
    alpha_diversity_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    primary_variable: str,
    output_file: str,
) -> pd.DataFrame:
    """Perform mannwhitneyu between every pair of groups for each alpha diversity metric.

    Similar to pairwise_alpha_diversity_calculations but without patient matching.
    Uses Mann-Whitney U test to compare all samples between groups directly.

    Args:
        alpha_diversity_df (pd.DataFrame): DataFrame containing alpha diversity metrics.
        metadata_df (pd.DataFrame): DataFrame containing sample metadata.
        primary_variable (str): The column in metadata_df to group by.
        output_file (str): Path to write out pandas dataframe as csv.

    Returns:
        pd.DataFrame: DataFrame containing pairwise comparison results with z-score scaled differences.
    """

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Initialize results list to store comparison results
    results_list = []

    # Merge alpha diversity with metadata
    merged_df = alpha_diversity_df.merge(
        metadata_df, left_on="sample_id", right_on="sample_id"
    )

    # Get groups from metadata
    groups = metadata_df[primary_variable].unique()

    for g1, g2 in combinations(groups, 2):
        comparison = f"{g1}_vs_{g2}"

        # Perform pairwise test for each alpha diversity metric
        for metric in [
            "Faith's phylogenetic diversity",
            "Observed Features",
            "Simpson's dominance index",
            "Simpson's diversity index",
        ]:
            group1_values = merged_df[merged_df[primary_variable] == g1][metric]
            group2_values = merged_df[merged_df[primary_variable] == g2][metric]

            # Skip if either group has no data
            if len(group1_values) == 0 or len(group2_values) == 0:
                print(f"No data for {metric} in groups {g1} or {g2}. Skipping.")
                continue

            # Perform Mann-Whitney U test
            stat, p_value = scipy.stats.mannwhitneyu(
                group1_values, group2_values, alternative="two-sided"
            )

            # Calculate z-score standardized difference
            pooled_values = pd.concat([group1_values, group2_values])
            pooled_std = pooled_values.std()
            raw_median_diff = group1_values.median() - group2_values.median()

            if pooled_std != 0:
                zscore_diff = raw_median_diff / pooled_std
            else:
                zscore_diff = 0.0  # Handle division by zero

            if raw_median_diff > 0:
                relative_change = "up in " + str(g1)
            elif raw_median_diff < 0:
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
                    "n_group1": len(group1_values),
                    "n_group2": len(group2_values),
                    "relative_change": relative_change,
                    "U-statistic": stat,
                    "p-value": p_value,
                    "median_diff": raw_median_diff,
                    "zscore_diff": zscore_diff,
                }
            )

    # Convert results list to DataFrame
    results_df = pd.DataFrame(results_list)

    # Helper function for significance stars
    def get_significance(p):
        if p < 0.0001:
            return "****"
        if p < 0.001:
            return "***"
        if p < 0.01:
            return "**"
        if p < 0.05:
            return "*"
        return "ns"

    # Add a column for significance based on p-value
    results_df["significance"] = results_df["p-value"].apply(get_significance)

    # Write results DataFrame to CSV
    results_df.to_csv(output_file, index=False)

    return results_df


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

    return ordination_results, distance_df
