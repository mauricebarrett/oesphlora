import os
import warnings
from typing import Optional

import biom  # type: ignore
import numpy as np
import pandas as pd
from biom.table import Table  # type: ignore
from gemelli.rpca import rpca  # type: ignore
from qiime2 import Artifact  # type: ignore
from qiime2.plugins.diversity.pipelines import alpha_phylogenetic
from qiime2.plugins.phylogeny.pipelines import align_to_tree_mafft_fasttree
from rpy2 import robjects as r
from rpy2.robjects import StrVector, pandas2ri, r
from skbio import DistanceMatrix, OrdinationResults
from skbio.diversity import alpha, beta_diversity
from sklearn.manifold import MDS


def generate_phylogenetic_tree_with_qiime2(input_fasta: str, threads, output_dir: str):
    """
    Generates a phylogenetic tree from an input FASTA file using QIIME2 MAFFT and FastTree.

    Parameters:
    - input_fasta (str): Path to the input FASTA file.
    - num_threads (int): Number of threads to use.
    - output_dir (str): Directory where output tree files should be saved.

    Returns:
    - qiime2.Artifact: The rooted tree QIIME2 artifact object.
    """

    rooted_tree_file = os.path.join(output_dir, "rooted_tree.qza")
    unrooted_tree_file = os.path.join(output_dir, "unrooted_tree.qza")
    # Skip if output directory already exists and read in trees
    if os.path.exists(output_dir):
        print(f"Output directory {output_dir} already exists.")
        print("Reading in existing trees...")
        # Load existing trees

        rooted_tree = Artifact.load(rooted_tree_file)
        unrooted_tree = Artifact.load(unrooted_tree_file)

        return rooted_tree, unrooted_tree
    print("Starting to generate phylogenetic tree...")

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Import sequences as a QIIME2 artifact
    rep_seq_artifact = Artifact.import_data(
        type="FeatureData[Sequence]", view=input_fasta
    )

    print("✅ Finished reading in sequences.")

    # Run alignment and tree-building
    action_results = align_to_tree_mafft_fasttree(
        sequences=rep_seq_artifact, n_threads=threads
    )

    # Extract rooted and unrooted tree objects
    rooted_tree = action_results.rooted_tree
    unrooted_tree = action_results.tree

    # Save trees to output directory

    rooted_tree.save(rooted_tree_file)
    unrooted_tree.save(unrooted_tree_file)
    print(f"Rooted tree saved to {rooted_tree_file}")
    print(f"Unrooted tree saved to {unrooted_tree_file}")

    return rooted_tree, unrooted_tree


def rarefy_df(
    df: pd.DataFrame,
    sampling_depth: Optional[int] = None,
    with_replacement: bool = False,
    random_seed: int = 1088,
) -> biom.Table:
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
    biom.Table
        The corresponding rarefied biom.Table object.
    """

    # Convert dataframe to biom Table
    biom_table = biom.Table(
        df.values,
        observation_ids=df.index.astype(str),
        sample_ids=df.columns.astype(str),
    )

    if sampling_depth is None:
        sampling_depth = int(min(biom_table.sum(axis="sample")))
    else:
        sampling_depth = int(sampling_depth)

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

    if rarefied_biom_table.is_empty():
        raise ValueError(
            "The rarefied table contains no samples or features. "
            "Verify your table is valid and that you provided a "
            "shallow enough sampling depth."
        )

        # Convert the biom table to a QIIME2 artifact
    rarefied_table_qza = Artifact.import_data(
        "FeatureTable[Frequency]", rarefied_biom_table
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

    return rarefied_df, rarefied_biom_table, rarefied_table_qza


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

    # skip if the output file already exists
    if os.path.exists(output_file):
        print(
            f"Output file {output_file} already exists. Skipping alpha diversity calculation."
        )
        return pd.read_csv(output_file)

    # Calculate alpha diversity metrics
    print("Calculating alpha diversity metrics...")
    alpha_diversity_metrics: dict[str, list] = {
        "sample_id": [],
        "Observed OTUs": [],
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
        alpha_diversity_metrics["Observed OTUs"].append(observed_otus)
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

    # Write out ordination and distance matrix to files
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    ordination.write(
        output_file.replace("distance_matrix.txt", "ordination.txt"),
        format="ordination",
    )
    robust_aitchison_dm.write(output_file)
    print(f"RPCA results saved to {output_file} and ordination file.")

    return ordination, robust_aitchison_dm


def perform_bray_curtis(
    rarefied_table_df: pd.DataFrame,
    output_file: str,
) -> DistanceMatrix:
    """Perform Bray-Curtis distance calculation on the rarefied OTU table.

    Args:
        rarefied_otu_table (pd.DataFrame): The rarefied OTU table where rows are taxa and columns are samples.

    Returns:
        DistanceMatrix: A Bray-Curtis distance matrix between the samples.
    """

    # skip if the output file already exists
    if os.path.exists(output_file):
        print(
            f"Output file {output_file} already exists. Skipping Bray-Curtis calculation."
        )
        # Return the existing DistanceMatrix
        return DistanceMatrix.read(output_file)

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

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    # Write the Bray-Curtis distance matrix
    bray_curtis_dm.write(output_file)
    print(f"Bray-Curtis distance matrix saved to {output_file}")

    return bray_curtis_dm


def perform_permanova_with_vegan(
    metadata: pd.DataFrame,
    distance_matrix: DistanceMatrix,
    primary_variable: str,
) -> tuple[dict, pd.DataFrame]:
    """Performs PERMANOVA analysis and tests homogeneity of dispersion using betadisper.

    Args:
        metadata (pd.DataFrame): A metadata table containing sample information.
        distance_matrix (DistanceMatrix): A distance matrix representing the distances between samples.

    Returns:
        tuple[dict, pd.DataFrame]: A dictionary containing PERMANOVA and betadisper results and the results DataFrame.
    """

    # Activate pandas-to-R conversion
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        pandas2ri.activate()

        # Convert pandas DataFrame to R dataframe
        r_metadata = pandas2ri.py2rpy(metadata)
        r_distance_matrix = pandas2ri.py2rpy(distance_matrix.to_data_frame())

    # Assign the metadata and distance matrix to the R environment
    r.assign("distance_matrix", r_distance_matrix)
    r.assign("metadata", r_metadata)
    r.assign("primary_variable", primary_variable)

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
