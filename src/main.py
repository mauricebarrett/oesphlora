import argparse
import glob
import io
import itertools
import os
import re
import subprocess
import sys
import uuid
import warnings
from collections import Counter
from pathlib import Path
from typing import List

import biom
import bionumpy as bnp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import qiime2
import seaborn as sns
from Bio import SeqIO
from biom.table import Table
from gemelli.rpca import rpca
from PIL import Image
from qiime2 import Artifact
from qiime2.plugins.diversity.pipelines import alpha_phylogenetic
from qiime2.plugins.feature_classifier.methods import classify_sklearn
from qiime2.plugins.phylogeny.pipelines import align_to_tree_mafft_fasttree
from rpy2.robjects import StrVector, pandas2ri, r
from rpy2.robjects.vectors import ListVector
from skbio import DistanceMatrix, OrdinationResults
from skbio.diversity import alpha, beta_diversity
from sklearn.manifold import MDS
from tqdm import tqdm


# Function to parse input arguments
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="A script to manage multiple language scripts with threading and directory support."
    )
    parser.add_argument(
        "-t", "--threads", type=int, required=True, help="Number of threads to use"
    )
    parser.add_argument(
        "-d", "--directory", type=str, required=True, help="Directory to work with"
    )
    args = parser.parse_args()
    return (
        args.threads,
        args.directory,
    )


# Function to gather QC data for a specific stage
def gather_qc_data(directory, stage, num_threads):
    # Print processing message
    print(
        f"Processing {stage} samples by:\n[1] Counting the number of reads per sample\n[2] Gleaning read length distribution per sample"
    )

    fastq_dir = os.path.join(directory, "fastq_files", stage)
    read_length_dir = os.path.join(directory, "qc", "read_length", stage)
    os.makedirs(read_length_dir, exist_ok=True)

    # Create a list to store read count data
    read_count_data = []

    # Check if there are any fastq files
    fastq_files = [
        f
        for f in os.listdir(fastq_dir)
        if f.endswith(".fastq") or f.endswith(".fastq.gz")
    ]
    # Sort fastq files to ensure consistency
    fastq_files.sort()

    # Check if there are any fastq files
    if fastq_files:
        print(f".fastq files found for {stage}, proceeding with processing.")
    else:
        print(f"No .fastq files found for {stage}. Exiting with error.")
        sys.exit(1)

    for filename in fastq_files:
        sample_name = os.path.splitext(filename)[0].split("_R")[0]
        print(f"Processing sample {sample_name}")

        # Counting number of reads in sample
        print(f"Counting number of reads in sample {sample_name}")
        fastq_file_path = os.path.join(fastq_dir, filename)
        count_cmd = ["seqkit", "stats", "-j", str(num_threads), fastq_file_path]
        count_cmd_output = subprocess.run(count_cmd, capture_output=True, text=True)
        count = int(count_cmd_output.stdout.splitlines()[1].split()[3].replace(",", ""))
        read_count_data.append({"sample_name": sample_name, "read_count": count})

        # Gleaning read length distribution in sample using Python
        if count > 50:
            print(f"Gleaning read length distribution in sample {sample_name}")
            read_length_file_path = os.path.join(
                read_length_dir, f"{sample_name}_read_length.txt"
            )

            # Run seqkit fx2tab -nl and capture output
            seqkit_cmd = [
                "seqkit",
                "fx2tab",
                "-j",
                str(num_threads),
                "-nl",
                fastq_file_path,
            ]
            seqkit_output = subprocess.run(
                seqkit_cmd, capture_output=True, text=True, check=True
            ).stdout

            # Extract and count read lengths
            read_lengths = [
                int(line.split("\t")[1]) for line in seqkit_output.splitlines()
            ]
            length_counts = Counter(read_lengths)

            # Write the counts to a file
            with open(read_length_file_path, "w") as f:
                for length, count in sorted(length_counts.items()):
                    f.write(f"{count}\t{length}\n")
        else:
            print(
                f"{sample_name} has insufficient reads ({count} reads)\nDeleting {stage} sample {sample_name}"
            )
            os.remove(fastq_file_path)

    # Save read count data as a CSV file
    if read_count_data:
        read_count_df = pd.DataFrame(read_count_data)
        read_count_csv_path = os.path.join(
            directory,
            "qc",
            "read_count",
            f"{stage}_reads_count.csv",
        )
        # Ensure the directory exists before saving
        os.makedirs(os.path.dirname(read_count_csv_path), exist_ok=True)
        # Save the CSV file
        read_count_df.to_csv(read_count_csv_path, index=False)


# Function to graph the quality scores of 16s FASTQ files
def graph_quality_scores(directory, stage):
    print(f"Graphing quality score of {stage} samples")

    # Get the path and read in files
    files_path = os.path.join(directory, "fastq_files", stage)
    fastq_files = [
        os.path.join(files_path, f)
        for f in os.listdir(files_path)
        if f.endswith(".fastq.gz")
    ]

    # Sort the fastq files to ensure consistency
    fastq_files.sort()

    # Prepare a figure to save all the quality plots
    output_dir = os.path.join(directory, "qc", "q_score", stage)
    os.makedirs(output_dir, exist_ok=True)

    # Start generating quality plots
    for fastq_file in fastq_files:
        # Extract sample name
        sample_name = os.path.basename(fastq_file).split(".")[0]

        print(f"Gathering Q score data for {sample_name}")

        # Define subset of fastq file to read
        subset_fastq_file = fastq_file.replace(".fastq.gz", "_subset.fastq")

        # Extract 2000 reads from the fastq file
        subprocess.run(
            f"seqkit head -n 1000 {fastq_file} > {subset_fastq_file}",
            shell=True,
            check=True,
        )

        # Read the fastq file
        f = bnp.open(subset_fastq_file)
        data = f.read()  # This reads the sequences and quality scores

        # delete the subset fastq file
        os.remove(subset_fastq_file)

        # Extract the quality scores
        quality_scores = data.quality

        # Determine the maximum read length
        read_length = max([len(read) for read in quality_scores])

        # Create an array to store quality scores for each position
        quality_matrix = np.full((len(quality_scores), read_length), np.nan)

        for i, read in enumerate(quality_scores):
            quality_matrix[i, : len(read)] = read

        # Calculate the statistics for each position
        mean_quality = np.nanmean(
            quality_matrix, axis=0
        )  # Mean quality score per position
        median_quality = np.nanmedian(
            quality_matrix, axis=0
        )  # Median quality score per position
        percentile_25 = np.nanpercentile(
            quality_matrix, 25, axis=0
        )  # 25th percentile per position
        percentile_75 = np.nanpercentile(
            quality_matrix, 75, axis=0
        )  # 75th percentile per position

        # Create the plot
        plt.figure(figsize=(14, 8))

        # Plotting median or mean with IQR
        plt.plot(
            range(1, read_length + 1),
            median_quality,
            label="Median Quality Score",
            color="blue",
            linewidth=2,
        )
        plt.fill_between(
            range(1, read_length + 1),
            percentile_25,
            percentile_75,
            color="skyblue",
            alpha=0.4,
            label="IQR (25th-75th Percentile)",
        )

        # Set quality score bands for visualization (optional)
        plt.axhspan(0, 20, color="red", alpha=0.1, label="Low Quality (<20)")
        plt.axhspan(20, 30, color="orange", alpha=0.1, label="Moderate Quality (20-30)")
        plt.axhspan(30, 50, color="green", alpha=0.1, label="High Quality (>30)")

        # Customize the plot
        plt.title(
            f"Quality Score Across Read Positions for Sample {sample_name}", fontsize=16
        )
        plt.xlabel("Base Position", fontsize=14)
        plt.ylabel("Phred Quality Score", fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.ylim(0, 45)
        plt.grid(visible=True, which="both", linestyle="--", linewidth=0.5)

        # Place the legend outside the plot
        plt.legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize=12)

        # Save each figure as a separate PDF file
        output_file = os.path.join(
            output_dir, f"{stage}_{sample_name}_phred_quality_lineplot.pdf"
        )
        plt.savefig(output_file, bbox_inches="tight")
        plt.close()

    print("Finished graphing")


# define funtion to remove primers using cutadapt


def remove_primers(
    directory,
    num_threads,
    adapter_fwd="CCTACGGGNGGCWGCAG",
    adapter_rev="GACTACHVGGGTATCTAATCC",
):
    """
    Removes primers from paired-end FASTQ files using cutadapt.

    Parameters:
        directory (str): Root directory containing fastq_files/demultiplexed/
        adapter_fwd (str): Forward adapter sequence
        adapter_rev (str): Reverse adapter sequence
    """
    files_path = os.path.join(directory, "fastq_files", "demultiplexed")
    output_path = os.path.join(directory, "fastq_files", "primers_removed")
    os.makedirs(output_path, exist_ok=True)

    r1_files = glob.glob(os.path.join(files_path, "*_R1_001.fastq.gz"))

    # Sort the fastq files to ensure consistency
    r1_files.sort()

    for r1 in r1_files:
        r2 = r1.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")

        if not os.path.exists(r2):
            print(f"⚠️  Missing R2 file for {os.path.basename(r1)}, skipping...")
            continue

        sample_name = os.path.basename(r1).replace("_R1_001.fastq.gz", "")
        out_r1 = os.path.join(output_path, f"{sample_name}_noprimers_R1.fastq.gz")
        out_r2 = os.path.join(output_path, f"{sample_name}_noprimers_R2.fastq.gz")

        cutadapt_cmd = [
            "cutadapt",
            "-j",
            str(num_threads),
            "--discard-untrimmed",
            "--overlap",
            "15",
            "--match-read-wildcards",
            "--minimum-length",
            "100",
            "-g",
            adapter_fwd,
            "-G",
            adapter_rev,
            "-o",
            out_r1,
            "-p",
            out_r2,
            r1,
            r2,
        ]

        print(f"🧪 Running cutadapt for sample: {sample_name}")
        result = subprocess.run(cutadapt_cmd, capture_output=True, text=True)

        if result.returncode == 0:
            print(f"✅ Done: {sample_name}")
        else:
            print(f"❌ Error with {sample_name}:\n{result.stderr}")


# Function to run the R script
def generate_asv_table(working_dir):
    try:
        subprocess.run(
            [
                "Rscript",
                "--working_dir",
                working_dir,
            ],
            check=True,
        )
        print("R script completed for 16S amplicon.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running R script for 16S: {e}")
        sys.exit(1)


def load_representative_sequences(rep_seq_fasta: str):
    """
    Load representative sequences from the specified FASTA file.

    Args:
        rep_seq_fasta (str): Path to the FASTA file containing representative sequences.

    Returns:
        Artifact: A QIIME2 Artifact object for the representative sequences.
    """

    print(f"Loading representative sequences from: {rep_seq_fasta}")

    # Load the sequence file as a QIIME2 Artifact
    rep_seq_artifact = Artifact.import_data(
        type="FeatureData[Sequence]", view=rep_seq_fasta
    )

    print("Successfully loaded representative sequences.")
    return rep_seq_artifact


def classify_taxonomy_nb(rep_seq_artifact, num_threads):
    """
    Perform taxonomic classification for 16S representative sequences.

    Args:
        rep_seq_artifact (Artifact): QIIME2 Artifact containing representative sequences.
        num_threads (int): Number of parallel jobs for classification.

    Returns:
        pd.DataFrame: DataFrame containing taxonomy classification results.
    """
    # Set the classifier path for 16S
    classifier_path = "/home/maurice/resources/marker_gene_databases/16s/greengenes/2024.09/v3_v4/qiime2/greengenes2_2024.09_v3v4_derep_classifier.qza"

    # Load the classifier artifact
    print(f"Loading classifier artifact from: {classifier_path}")
    classifier_artifact = Artifact.load(classifier_path)

    # Perform taxonomic classification
    print("Classifying 16S representative sequences...")
    taxonomy_classification = classify_sklearn(
        classifier=classifier_artifact,
        reads=rep_seq_artifact,
        confidence=0.7,
        n_jobs=2,
    )

    print("Taxonomic classification completed.")

    # Convert taxonomy to a pandas DataFrame
    taxonomy_df = taxonomy_classification.classification.view(pd.DataFrame)
    print("Taxonomic classification completed successfully.")

    return taxonomy_df


def reformat_taxonomy(taxonomy_df: pd.DataFrame):
    """
    Reformat the taxonomy DataFrame to split taxonomic ranks into separate columns.
    Args:
        taxonomy_df (pd.DataFrame): DataFrame containing taxonomy classification results.
    Returns:
        pd.DataFrame: Reformatted taxonomy DataFrame with separate columns for each taxonomic rank.
    """
    # Split 'Taxon' column into separate columns
    taxonomy_df[
        ["domain", "phylum", "class", "order", "family", "genus", "species"]
    ] = taxonomy_df["Taxon"].str.split(";", expand=True)

    # Clean data naming
    columns_to_clean = [
        "domain",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]
    for col in columns_to_clean:
        taxonomy_df[col] = taxonomy_df[col].str.strip()
        taxonomy_df[col] = taxonomy_df[col].replace(
            r"^d__$|^k__$|^p__$|^c__$|^o__$|^f__$|^g__$|^s__$",
            "unclassified",
            regex=True,
        )
        taxonomy_df[col] = taxonomy_df[col].fillna(value="unclassified")
        taxonomy_df[col] = taxonomy_df[col].replace(
            to_replace="None", value="unclassified"
        )
        taxonomy_df[col] = taxonomy_df[col].replace(
            r"^d__|^k__|^p__|^c__|^o__|^f__|^g__|^s__", "", regex=True
        )
        taxonomy_df[col] = taxonomy_df[col].replace(
            to_replace="Unassigned", value="unclassified"
        )

    taxonomy_df["Taxon"] = taxonomy_df["Taxon"].str.replace("; ", ";", regex=False)
    taxonomy_df["Taxon"] = taxonomy_df["Taxon"].replace(
        r";p__|;c__|;o__|;f__|;g__|;s__", ";", regex=True
    )
    taxonomy_df["Taxon"] = taxonomy_df["Taxon"].replace(r"d__|k__", "", regex=True)
    taxonomy_df["Taxon"] = taxonomy_df["Taxon"].replace(
        to_replace="Unassigned", value="unclassified"
    )

    return taxonomy_df


def filter_taxonomy_features(taxonomy_df, keywords: str):
    """
    Filter features from taxonomy DataFrame based on keywords and save the filtered data.

    Args:
        taxonomy_df (pd.DataFrame): DataFrame containing taxonomy data.
        keywords (list): List of keywords to filter out (case-insensitive).

    Returns:
        tuple: Filtered taxonomy DataFrame and list of feature IDs to filter out.
    """
    # Create a regex pattern from keywords
    pattern = "|".join(keywords)

    # Identify rows containing specified keywords
    features_to_filter = taxonomy_df[
        taxonomy_df.astype(str)
        .apply(lambda col: col.str.contains(pattern, case=False, na=False))
        .any(axis=1)
    ]

    # Extract 'Feature ID' values to filter out
    feature_ids_to_filter = features_to_filter.index.tolist()

    # Filter taxonomy_df to exclude these features
    filtered_taxonomy_df = taxonomy_df.drop(index=feature_ids_to_filter)

    return filtered_taxonomy_df, feature_ids_to_filter


def filter_asv_table(asv_table_df: pd.DataFrame, feature_ids_to_filter: list):
    """
    Read and filter the ASV table based on feature IDs to exclude.

    Args:
        asv_table_df (pd.DataFrame): ASV table DataFrame.
        feature_ids_to_filter (list): List of feature IDs to remove.

    Returns:
        pd.DataFrame: Filtered ASV table.
    """

    # Filter ASV table to exclude features identified for removal
    filtered_asv_table_df = asv_table_df.drop(
        index=feature_ids_to_filter, errors="ignore"
    )

    # Calculate the total count for both dataframes
    original_count = asv_table_df.sum().sum()
    filtered_count = filtered_asv_table_df.sum().sum()

    # Calculate percentage of counts lost
    counts_lost = original_count - filtered_count
    percentage_lost = (counts_lost / original_count) * 100
    print(f"Percentage of count lost: {percentage_lost:.2f}%")

    return filtered_asv_table_df


def filter_fasta_sequences(wor_dir, feature_ids_to_filter):
    """
    Filter sequences from a FASTA file based on feature IDs to exclude.

    Args:
        wor_dir (str): The working directory containing the representative sequences.
        feature_ids_to_filter (list): List of feature IDs to remove.
    """
    # Construct input and output FASTA paths
    input_fasta = os.path.join(
        wor_dir,
        "representative_sequences",
        "asv_sequences.fasta",
    )
    output_fasta = os.path.join(
        wor_dir,
        "representative_sequences",
        "asv_filtered_sequences.fasta",
    )

    # Ensure the input FASTA file exists
    if not Path(input_fasta).exists():
        raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")

    # Create the output directory if it doesn't exist
    Path(os.path.dirname(output_fasta)).mkdir(parents=True, exist_ok=True)

    # Open the output FASTA file and filter sequences
    with open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id not in feature_ids_to_filter:
                SeqIO.write(record, output_handle, "fasta")

    print(f"Filtered FASTA saved to {output_fasta}")


def generate_phylogenetic_tree_with_qiime2(input_fasta: str, num_threads):
    """
    Generates a phylogenetic tree from an input FASTA file using QIIME2 MAFFT and FastTree.

    Parameters:
    - input_fasta (str): Path to the input FASTA file.
    - output_dir (str): Directory where output tree files should be saved.

    Returns:
    - qiime2.Artifact: The rooted tree QIIME2 artifact object.
    """

    # Import sequences as a QIIME2 artifact
    rep_seq_artifact = qiime2.Artifact.import_data(
        type="FeatureData[Sequence]", view=input_fasta
    )

    print("✅ Finished reading in sequences.")

    # Run alignment and tree-building
    action_results = align_to_tree_mafft_fasttree(
        sequences=rep_seq_artifact, n_threads=num_threads
    )

    # Extract rooted tree object
    rooted_tree = action_results.rooted_tree
    unrooted_tree = action_results.tree

    return rooted_tree, unrooted_tree


def create_taxa_tables(
    filtered_asv_table_df: pd.DataFrame,
    filtered_taxonomy_df: pd.DataFrame,
    rank: str,
):
    """
    Create taxa tables for species, genus, and family ranks and save them as CSV files.

    Args:
        filtered_asv_table_df (pd.DataFrame): Feature Filtered ASV table DataFrame.
        filtered_taxonomy_df (pd.DataFrame): Feature Filtered taxonomy DataFrame.
        rank (str): The taxonomic rank to process ('species', 'genus', or 'family').
    """

    # Initialize dictionary to store grouped data
    taxa_table_dict = {}

    # Group and sum ASV abundances by taxonomic rank
    for col in tqdm(
        filtered_asv_table_df.columns, desc=f"Processing columns for {rank}"
    ):
        grouped_sums = (
            filtered_asv_table_df[col].groupby(filtered_taxonomy_df[rank]).sum()
        )
        taxa_table_dict[col] = grouped_sums

    # Convert to DataFrame
    return pd.DataFrame(taxa_table_dict)


def depth_filtering(table: pd.DataFrame, depth_threshold: float) -> pd.DataFrame:
    """Remove samples from feature table with read depth below threshold."""
    sample_depths = table.sum(axis=0)

    # Remove samples with read depth below threshold
    return table.loc[:, sample_depths >= depth_threshold]


def rarefy_df(
    df: pd.DataFrame,
    with_replacement: bool = False,
    random_seed: int = 1088,
) -> pd.DataFrame:
    """
    Rarefy a pandas DataFrame after converting it to biom.Table and then back.

    Parameters:
    ----------
    df : pd.DataFrame
        Input dataframe with samples as columns and features as rows.

    sampling_depth : int, optional
        Rarefaction depth. If None, automatically set to the minimum column sum.

    with_replacement : bool, default False
        Perform rarefaction sampling with or without replacement.

    random_seed : int, default 1088
        Seed for random number generation to ensure reproducibility.

    Returns:
    -------
    pd.DataFrame
        Rarefied dataframe in the same format as input.
    """

    # Calculate minimum sample depth
    sampling_depth = df.sum(axis=0).min()

    # Convert dataframe to biom Table
    biom_table = biom.Table(
        df.values,
        observation_ids=df.index.astype(str),
        sample_ids=df.columns.astype(str),
    )

    # Filter out samples that don't meet sampling_depth if sampling with replacement
    biom_table = biom_table.filter(
        lambda v, i, m: v.sum() >= sampling_depth, inplace=False, axis="sample"
    )

    # Rarefy biom Table
    rarefied_biom_table = biom_table.subsample(
        sampling_depth,
        axis="sample",
        by_id=False,
        with_replacement=with_replacement,
        seed=random_seed,
    )

    # Check if rarefied table is empty
    if rarefied_biom_table.is_empty():
        raise ValueError(
            "The rarefied table contains no samples or features. "
            "Verify your table is valid and that you provided a "
            "shallow enough sampling depth."
        )

    # Convert biom Table back to pandas DataFrame
    rarefied_df = pd.DataFrame(
        rarefied_biom_table.matrix_data.toarray(),
        index=rarefied_biom_table.ids(axis="observation"),
        columns=rarefied_biom_table.ids(axis="sample"),
    )

    return rarefied_df, rarefied_biom_table


def calculate_alpha_diversity(normalized_table: pd.DataFrame) -> pd.DataFrame:
    """Calculates alpha diversity metrics for each sample in the normalized table.

    Args:
        normalized_table (pd.DataFrame): The normalized OTU table where rows are taxa and columns are samples.

    Returns:
        pd.DataFrame: A DataFrame containing alpha diversity metrics for each sample.
    """
    alpha_diversity_metrics: dict[str, list] = {
        "sample_id": [],
        "Observed OTUs": [],
        "Simpson's dominance index": [],
        "Simpson's diversity index": [],
    }

    for sample in normalized_table.columns:
        # Calculate observed OTUs
        observed_otus = alpha.observed_otus(normalized_table[sample])

        # Calculate Simpson's dominance index
        dominance = alpha.dominance(normalized_table[sample])

        # Calculate Simpson diversity index
        simpson = alpha.simpson(normalized_table[sample])

        alpha_diversity_metrics["sample_id"].append(sample)
        alpha_diversity_metrics["Observed OTUs"].append(observed_otus)
        alpha_diversity_metrics["Simpson's dominance index"].append(dominance)
        alpha_diversity_metrics["Simpson's diversity index"].append(simpson)

    alpha_diversity_df = pd.DataFrame(alpha_diversity_metrics)

    return alpha_diversity_df


def calculate_phylo_alpha_diversity_qiime2(
    normalized_table: pd.DataFrame, rooted_tree, metric: str = "faith_pd"
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


def prevalance_filtering(
    otu_table: pd.DataFrame, prevalence_threshold: float
) -> pd.DataFrame:
    """Remove features from feature table with prevalence below threshold."""
    feature_prevalence = (otu_table > 0).sum(axis=1) / otu_table.shape[1] * 100

    # Remove features with prevalence below threshold
    return otu_table.loc[feature_prevalence >= prevalence_threshold, :]


def linda_daa(
    asv_table: pd.DataFrame,
    metadata: pd.DataFrame,
    primary_variable: str,
    taxonomy_table: pd.DataFrame,
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

        # Check if the primary variable is categorical
        metadata_subset["Diagnosis"] = metadata_subset[
            "Diagnosis"
        ].cat.remove_unused_categories()

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

        if not results:
            print("No results to concatenate. All comparisons were skipped.")
            return pd.DataFrame(), pd.DataFrame()  # Return empty DataFrames

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

        daa_results_taxa_fil_df = daa_results_taxa_df[
            daa_results_taxa_df["padj"] < 0.05
        ]

    return daa_results_taxa_df, daa_results_taxa_fil_df


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
    r_script_path = os.path.join(current_script_dir, "R/alpha_diversity_boxplots.R")
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


def perform_bray_curtis(
    rarefied_otu_table: pd.DataFrame,
) -> DistanceMatrix:
    """Perform Bray-Curtis distance calculation on the rarefied OTU table.

    Args:
        rarefied_otu_table (pd.DataFrame): The rarefied OTU table where rows are taxa and columns are samples.

    Returns:
        DistanceMatrix: A Bray-Curtis distance matrix between the samples.
    """
    # Transpose the OTU table to have samples as rows and taxa as columns (as required by Bray-Curtis)
    otu_data = rarefied_otu_table.T.to_numpy()

    # Perform Bray-Curtis calculation
    bray_curtis_results = beta_diversity(
        metric="braycurtis", counts=otu_data, ids=rarefied_otu_table.columns
    )

    # Convert to a DistanceMatrix object for further downstream analysis
    bray_curtis_dm = DistanceMatrix(
        bray_curtis_results.data, ids=bray_curtis_results.ids
    )

    # Convert to data-frame
    bray_curtis_df = bray_curtis_dm.to_data_frame()

    return bray_curtis_dm, bray_curtis_df


def perform_permanova(
    metadata: pd.DataFrame, distance_matrix: DistanceMatrix
) -> tuple[dict, pd.DataFrame]:
    """Performs PERMANOVA analysis and tests homogeneity of dispersion using betadisper.

    Args:
        metadata (pd.DataFrame): A metadata table containing sample information.
        distance_matrix (DistanceMatrix): A distance matrix representing the distances between samples.

    Returns:
        tuple[dict, pd.DataFrame]: A dictionary containing PERMANOVA and betadisper results and the results DataFrame.
    """
    # Convert UUID columns to string if any
    for column in metadata.columns:
        if metadata[column].iloc[0].__class__ == uuid.UUID:
            metadata[column] = metadata[column].apply(str)

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

    # Run PERMANOVA and betadisper using R code
    r("""
    library(vegan)
      
    distance_matrix<-as.dist(distance_matrix)
      
    # Perform PERMANOVA
    adonis_result <- adonis2(distance_matrix ~internal_treatment_code, data = metadata, permutations = 999)

    # Extract specific PERMANOVA results
    pseudo_F <- adonis_result$`F`[1]
    p_value <- adonis_result$`Pr(>F)`[1]
    r_squared <- adonis_result$R2[1]

    # Perform betadisper to assess homogeneity of dispersion
    dispersion_result <- betadisper(distance_matrix, metadata$internal_treatment_code)
    dispersion_test <- permutest(dispersion_result, permutations = 999)

    # Extract betadisper results
    dispersion_F <- dispersion_test$tab$`F`[1]
    dispersion_p_value <- dispersion_test$tab$`Pr(>F)`[1]
      
    results_df <- data.frame(
    metric = c("pseudo_F", "p_value", "r_squared", "dispersion_F", "dispersion_p_value"),
    value = c(pseudo_F, p_value, r_squared, dispersion_F, dispersion_p_value)
    )
    """)

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


def plot_nmds_with_ggplot2(
    nmds_coordinates: pd.DataFrame,
    metadata: pd.DataFrame,
    primary_variable: str,
    stress_value: float,
    permanova_results: dict[str, float],
) -> None:
    """
    Save NMDS coordinates and metadata, then call the R script.

    Args:
        nmds_coordinates (pd.DataFrame): DataFrame containing NMDS coordinates with samples as rows.
        metadata (pd.DataFrame): Metadata for samples.
        primary_variable (str): Column in metadata used for coloring the plot.
        stress_value (float): NMDS stress value indicating the goodness of fit.
        permanova_results (dict[str, float]): PERMANOVA results
    """

    # Convert UUID columns to strings
    uuid_columns = [
        "sample_set_id",
        "sequencing_run_id",
        "location_id",
        "site_id",
        "customer_id",
    ]
    for col in uuid_columns:
        if col in metadata.columns:
            metadata[col] = metadata[col].astype(str)

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
    r_script_path = os.path.join(current_script_dir, "../r_scripts/nmds_plot.R")
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


def perform_rpca(
    biom_table: Table,
) -> tuple[OrdinationResults, DistanceMatrix, pd.DataFrame]:
    """Perform Robust Aitchison PCA (RPCA) on BIOM table.

    Args:
        table (Table): BIOM table

    Returns:
        tuple[OrdinationResults, DistanceMatrix]: coordinates and distance matrix
    """
    # Preform RPCA
    ordination, robust_aitchison_distmat = rpca(biom_table, min_sample_count=500)

    # Convert to Pandas DataFrame
    robust_aitchison_df = robust_aitchison_distmat.to_data_frame()

    return ordination, robust_aitchison_distmat, robust_aitchison_df


def copy_to_dropbox(source: str, destination: str, exclude_dirs: List[str]):
    """
    Copy files from source to destination using rclone, excluding specific directories.

    Args:
        source (str): Source directory to copy from.
        destination (str): Destination directory to copy to.
        exclude_dirs (List[str]): List of directories to exclude from the copy operation.
    """

    # Base rclone command
    cmd = [
        "rclone",
        "copy",
        source,
        destination,
        "--tpslimit",
        "10",
        "--retries",
        "4",
        "--stats-one-line",
        "--progress",
    ]

    # Add exclude flags
    for exclude in exclude_dirs:
        cmd.extend(["--exclude", f"{exclude}/**"])

    print(f"Copying files from {source} to {destination}...")
    print(f"Excluding directories: {', '.join(exclude_dirs)}")
    # Run the command
    result = subprocess.run(cmd)

    if result.returncode != 0:
        print("Dropbox copy failed! Exiting...")
        sys.exit(1)


def main():
    # Parse the arguments using the parse_arguments function
    num_threads, directory = parse_arguments()

    # Keywords to filter
    keywords_to_filter = ["Chloroplast", "Mitochondrion"]

    # Display the values for now (can be used later in the script)
    print(f"Number of threads: {num_threads}")
    print(f"Working directory: {directory}")

    # # Gather QC data for each stage
    # gather_qc_data(directory, "demultiplexed", num_threads)

    # # Graph quality scores of 16s FASTQ files
    # graph_quality_scores(directory, "demultiplexed")

    # Remove primers from demultiplexed FASTQ files
    # remove_primers(directory, num_threads)

    # Gather QC data for each stage
    # gather_qc_data(directory, "primers_removed", num_threads)

    # # Graph quality scores of 16s FASTQ files
    # graph_quality_scores(directory, "primers_removed")

    # # Run the R script to generate ASV table
    # generate_asv_table(directory)

    # # Define the path to the representative sequences fasta file
    # rep_seq_fasta = os.path.join(
    #     directory,
    #     "representative_sequences",
    #     "asv_sequences.fasta",
    # )

    # # Load representative sequences
    # rep_seq_artifact = load_representative_sequences(rep_seq_fasta)

    # # Use a pre-trained classifier
    # taxonomy_df = classify_taxonomy_nb(rep_seq_artifact, num_threads)

    # # Reformat taxonomy
    # taxonomy_df = reformat_taxonomy(taxonomy_df)

    # # Create directory for taxonomic classification output
    tc_dir = os.path.join(directory, "taxonomic_classification")
    # Path(tc_dir).mkdir(parents=True, exist_ok=True)

    # # Create path for initial taxonomic classification
    # tc_file_path = os.path.join(tc_dir, "taxonomy_data.csv")

    # # Write out initial taxonomic classification file
    # taxonomy_df.to_csv(tc_file_path, index=True)

    # # Filter features from taxonomy
    # filtered_taxonomy_df, feature_ids_to_filter = filter_taxonomy_features(
    #     taxonomy_df, keywords_to_filter
    # )

    # # Save filtered taxonomy DataFrame
    filtered_taxonomy_file = os.path.join(tc_dir, "taxonomy_data_filtered.csv")
    # filtered_taxonomy_df.to_csv(filtered_taxonomy_file, index=True)
    # print(f"Filtered taxonomy saved to: {filtered_taxonomy_file}")

    # # Path to ASV table
    asv_table_dir = os.path.join(directory, "tables", "asv_tables")
    # asv_table_file = os.path.join(asv_table_dir, "asv_table.csv")

    # # Read ASV table
    # print(f"Reading ASV table from: {asv_table_file}")
    # asv_table_df = pd.read_csv(asv_table_file, index_col=0)

    # # Filter ASV table
    # asv_table_feat_fil_df = filter_asv_table(asv_table_df, feature_ids_to_filter)

    # # Save filtered ASV table
    asv_table_feat_fil_file = os.path.join(
        asv_table_dir, "asv_table_feature_filtered.csv"
    )
    # asv_table_feat_fil_df.to_csv(asv_table_feat_fil_file, index=True)
    # print(f"Filtered ASV table saved to: {asv_table_feat_fil_file}")

    # # Filter FASTA sequences
    # filter_fasta_sequences(directory, feature_ids_to_filter)

    # Define directory for phlogenetics anaylsis
    phylo_dir = os.path.join(directory, "phylogenetic_tree")

    rooted_tree_qza_file = os.path.join(phylo_dir, "rooted_tree.qza")

    if not os.path.exists(rooted_tree_qza_file):
        # Generate phylogenetic tree using QIIME2
        rooted_tree, unrooted_tree = generate_phylogenetic_tree_with_qiime2(
            os.path.join(
                directory, "representative_sequences", "asv_filtered_sequences.fasta"
            ),
            num_threads,
        )

        Path(phylo_dir).mkdir(parents=True, exist_ok=True)

        # Save the rooted tree
        rooted_tree_file = os.path.join(phylo_dir, "rooted_tree.qza")
        rooted_tree.save(rooted_tree_file)
        # Save the unrooted tree
        unrooted_tree_file = os.path.join(phylo_dir, "unrooted_tree.qza")
        unrooted_tree.save(unrooted_tree_file)
    else:
        print(f"Rooted tree already exists at {rooted_tree_qza_file}")
        # Load the rooted tree
        rooted_tree = Artifact.load(rooted_tree_qza_file)
        print("Rooted tree loaded successfully.")

    # Read in the filtered ASV table
    asv_table_feat_fil_df = pd.read_csv(asv_table_feat_fil_file, index_col=0)

    # Remove "OAC182-1", "OAC182-2","OAC182-3","OAC182-4","OAC182-5"
    asv_table_feat_fil_df = asv_table_feat_fil_df.drop(
        ["OAC182-1", "OAC182-2", "OAC182-3", "OAC182-4", "OAC182-5"], axis=1
    )

    # Remove saples that dont start with "OAC"
    asv_table_feat_fil_df = asv_table_feat_fil_df.filter(regex="^OAC", axis=1)
    # Remove samples that start with "OAC" but not "OAC182"

    # Filter OTU table based on read depth
    asv_table_dep_fil_df = depth_filtering(asv_table_feat_fil_df, depth_threshold=2000)

    # Read in the filtered taxonomy table
    filtered_taxonomy_df = pd.read_csv(filtered_taxonomy_file, index_col=0)

    # Directory for taxa tables
    taxa_tables_dir = os.path.join(directory, "tables", "taxa_tables")
    Path(taxa_tables_dir).mkdir(parents=True, exist_ok=True)

    # Taxonomic ranks to process
    ranks = ["species", "genus", "family"]

    for rank in ranks:
        print(f"Processing taxonomic rank: {rank}")

        # Create taxa tables
        taxa_table_df = create_taxa_tables(
            asv_table_feat_fil_df,
            filtered_taxonomy_df,
            rank=rank,
        )

        # Save taxa table to CSV
        output_file = os.path.join(taxa_tables_dir, f"{rank}_table.csv")
        taxa_table_df.to_csv(output_file, index=True)
        print(f"Taxa table for {rank} saved to: {output_file}")

    # Define path to metadata file
    metadata_file = os.path.join(directory, "metadata", "metadata.csv")

    # Read metadata file
    metadata_df = pd.read_csv(metadata_file)
    # Rename the first column to 'patient_id'
    metadata_df.rename(columns={metadata_df.columns[0]: "patient_id"}, inplace=True)

    diagnosis_map = {
        "1": "Healthy",
        "2": "GORD",
        "3": "Barrett's Esophagus",
        "4": "Dysplasia",
        "5": "Esophageal adenocarcinoma",
        "6": "Metastatic Esophageal adenocarcinoma",
    }

    metadata_df["Diagnosis"] = metadata_df["Diagnosis"].replace(diagnosis_map)

    # Extract sample names from the ASV table as a dataframe
    sample_data = asv_table_dep_fil_df.columns.to_frame(name="sample_name")
    # Get patient IDs from sample names by removing evething after the first -
    sample_data["patient_id"] = sample_data["sample_name"].str.extract(r"^([^-]+)")
    # Get sample locations from sample names by removing everything before the first -
    sample_data["sample_location"] = sample_data["sample_name"].str.extract(r"-(.*)")

    # Merge sample data with metadata
    merged_df = pd.merge(
        sample_data,
        metadata_df,
        left_on="patient_id",
        right_on="patient_id",
        how="inner",
    )

    # Set sample_id as the index
    merged_df.set_index("sample_name", inplace=True)

    # Reorder the levels of the Diagnosis column
    ordered_levels = [
        "Healthy",
        "GORD",
        "Barrett's Esophagus",
        "Dysplasia",
        "Esophageal adenocarcinoma",
        "Metastatic Esophageal adenocarcinoma",
    ]

    merged_df["Diagnosis"] = pd.Categorical(
        merged_df["Diagnosis"], categories=ordered_levels, ordered=True
    )
    # Sort the DataFrame by the Diagnosis column
    merged_df.sort_values(by="Diagnosis", inplace=True)

    # Save out merged data as csv
    merged_data_file = os.path.join(directory, "metadata", "merged_data.csv")

    # Save the merged DataFrame to a CSV file
    merged_df.to_csv(merged_data_file, index=True)

    ##########################################
    # Caculate Alpha diversity all samples
    ##########################################

    print("Performing Alpha-diversity analysis")

    # rarify the ASV table
    asv_table_rarefied, rarefied_biom_table = rarefy_df(
        asv_table_dep_fil_df,
        with_replacement=False,
        random_seed=1088,
    )

    # Calculate alpha diversity metrics
    alpha_diversity_df = calculate_alpha_diversity(asv_table_rarefied)

    # Convert the biom table to a QIIME2 artifact
    rarefied_table_qza = Artifact.import_data(
        "FeatureTable[Frequency]", rarefied_biom_table
    )

    # Read in rootee tree as QIIME2 artifact
    rooted_tree_qza = Artifact.load(rooted_tree_qza_file)

    # Calculate phylogenetic alpha diversity
    phylo_alpha_metrics_df = calculate_phylo_alpha_diversity_qiime2(
        rarefied_table_qza, rooted_tree_qza
    )

    alpha_diversity_merged_df = phylo_alpha_metrics_df.merge(
        alpha_diversity_df, left_on="sample_id", right_on="sample_id"
    )

    # Define diversity metircs folder
    div_met_dir = os.path.join(directory, "diversity_metrics")

    # Define dir to output alpha results
    alpha_div_dir = os.path.join(div_met_dir, "alpha_diversity")

    # Create the directory if it doesn't exist
    os.makedirs(alpha_div_dir, exist_ok=True)

    # Create file path to save the results
    alpha_diversity_file_path = os.path.join(
        alpha_div_dir,
        "alpha_diversity_metrics.csv",
    )

    # Save the alpha diversity metrics to a CSV file
    alpha_diversity_merged_df.to_csv(alpha_diversity_file_path, index=False)

    print(f"Alpha diversity metrics saved to: {alpha_diversity_file_path}")

    #######################################################
    #
    #######################################################

    locations = merged_df["sample_location"].unique()

    for location in locations:
        print(f"Processing location: {location}")

        # Filter the merged DataFrame for the current location
        location_df = merged_df[merged_df["sample_location"] == location]

        # Filter the ASV table for the current location
        asv_table_location = asv_table_dep_fil_df[location_df.index]

        # Check if there are enough samples to perform analysis
        if asv_table_location.shape[1] < 4:
            print(
                f"Skipping {location}: only {asv_table_location.shape[1]} samples found (minimum 4 required)"
            )
            continue

        # Run LINDA analysis
        daa_results_taxa_df, daa_results_taxa_fil_df = linda_daa(
            asv_table_location,
            location_df,
            primary_variable="Diagnosis",
            taxonomy_table=filtered_taxonomy_df,
        )

        # Define directory for DAA results
        daa_dir = os.path.join(directory, "daa_results", location)
        # Create the directory if it doesn't exist
        os.makedirs(daa_dir, exist_ok=True)

        # Define path to save the full results
        full_results_file = os.path.join(
            daa_dir, f"diagnosis_daa_results_location_{location}.csv"
        )

        # Save the full results
        daa_results_taxa_df.to_csv(full_results_file, index=False)
        print(f"Full results saved to: {full_results_file}")

        # Define path to save the filtered results
        filtered_results_file = os.path.join(
            daa_dir, f"diagnosis_daa_results_filtered_location_{location}.csv"
        )
        # Save the filtered results
        daa_results_taxa_fil_df.to_csv(filtered_results_file, index=False)
        print(f"Filtered results saved to: {filtered_results_file}")

        ######################
        # Alpha Diversity
        ######################

        # Subset the alpha diversity DataFrame for the current location
        alpha_diversity_location_df = alpha_diversity_merged_df[
            alpha_diversity_merged_df["sample_id"].isin(location_df.index)
        ]

        # Plot alpha diversity boxplots for the current location
        img = plot_alpha_diversity_boxplots_with_ggplot2(
            alpha_diversity_location_df, location_df, "Diagnosis"
        )

        # Define directory for figures
        alpha_figures_dir = os.path.join(directory, "figures", "alpha_diversity")
        # Create the directory if it doesn't exist
        os.makedirs(alpha_figures_dir, exist_ok=True)

        # Save the plot as a PDF
        plot_file_path = os.path.join(
            alpha_figures_dir, f"alpha_diversity_boxplots_location_{location}.png"
        )

        # Save the image as png
        img.save(plot_file_path, format="PNG")
        print(f"Alpha diversity boxplot saved to: {plot_file_path}")

        #########################
        # Beta Diversity
        #########################

        # define directory for beta diversity
        beta_div_dir = os.path.join(directory, "diversity_metrics", "beta_diversity")

        # Create the directory if it doesn't exist
        os.makedirs(beta_div_dir, exist_ok=True)

        # rarify the ASV table
        rarefied_asv_table_location, rarefied_biom_table = rarefy_df(
            asv_table_location,
            with_replacement=False,
            random_seed=1088,
        )

        # Perform Bray-Curtis distance calculation
        bray_curtis_dm, bray_curtis_df = perform_bray_curtis(
            rarefied_asv_table_location
        )

        # Save the Bray-Curtis distance matrix to a CSV file
        bray_curtis_file_path = os.path.join(
            beta_div_dir, f"bray_curtis_distance_matrix_location_{location}.csv"
        )

        bray_curtis_df.to_csv(bray_curtis_file_path, index=False)

        print(f"Bray-Curtis distance matrix saved to: {bray_curtis_file_path}")

        # Perform RPCA
        ordination, robust_aitchison_dm, robust_aitchison_df = perform_rpca(
            rarefied_biom_table
        )

        # Save the RPCA results to a CSV file
        rpca_file_path = os.path.join(
            beta_div_dir, f"rpca_distance_matrix_location_{location}.csv"
        )
        robust_aitchison_df.to_csv(rpca_file_path, index=False)
        print(f"RPCA distance matrix saved to: {rpca_file_path}")

    # Copy the results to Dropbox
    copy_to_dropbox(
        source=directory,
        destination="onedrive:/sharing/oesphlora/bioinformatics",
        exclude_dirs=["data", "fastq_files"],
    )

    # # Run LINDA analysis
    # daa_results_taxa_df, daa_results_taxa_fil_df = linda_daa(
    #     asv_table_dep_fil_df,
    #     merged_df,
    #     primary_variable="Diagnosis",
    #     taxonomy_table=filtered_taxonomy_df,
    # )

    # # Define directory for DAA results
    # daa_dir = os.path.join(directory, "tables", "daa_results")
    # # Create the directory if it doesn't exist
    # os.makedirs(daa_dir, exist_ok=True)

    # # Define path to save the full results
    # full_results_file = os.path.join(daa_dir, "full_results.csv")

    # # Save the full results
    # daa_results_taxa_df.to_csv(full_results_file, index=False)
    # print(f"Full results saved to: {full_results_file}")

    # # Define path to save the filtered results
    # filtered_results_file = os.path.join(daa_dir, "filtered_results.csv")
    # # Save the filtered results
    # daa_results_taxa_fil_df.to_csv(filtered_results_file, index=False)
    # print(f"Full results saved to: {full_results_file}")


if __name__ == "__main__":
    main()
