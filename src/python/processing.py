import glob
import os
import subprocess
import sys
from typing import Optional

import pandas as pd
from qiime2 import Artifact
from qiime2.plugins import feature_classifier


def remove_primers(
    input_directory,
    num_threads,
    output_directory,
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

    os.makedirs(output_directory, exist_ok=True)

    r1_files = glob.glob(os.path.join(input_directory, "*_R1_001.fastq.gz"))

    # Sort the fastq files to ensure consistency
    r1_files.sort()

    for r1 in r1_files:
        r2 = r1.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")

        if not os.path.exists(r2):
            print(f"⚠️  Missing R2 file for {os.path.basename(r1)}, skipping...")
            continue

        sample_name = os.path.basename(r1).replace("_R1_001.fastq.gz", "")
        out_r1 = os.path.join(output_directory, f"{sample_name}_noprimers_R1.fastq.gz")
        out_r2 = os.path.join(output_directory, f"{sample_name}_noprimers_R2.fastq.gz")

        # skip if output files already exist
        if os.path.exists(out_r1) and os.path.exists(out_r2):
            print(f"✅ Output files already exist for {sample_name}, skipping...")
            continue

        cutadapt_cmd = [
            "cutadapt",
            "-j",
            str(num_threads),
            "--discard-untrimmed",
            "--overlap",
            "12",
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
def generate_asv_table(wor_dir):
    # Define path to asv table
    asv_table_path = os.path.join(wor_dir, "tables", "asv_tables")

    # Check if the asv table exits, if it does, skip the R script
    if os.path.exists(asv_table_path):
        print(f"✅ ASV table already exists at {asv_table_path}, skipping R script.")
        return

    # Run the R script to generate ASV table
    try:
        subprocess.run(
            [
                "Rscript",
                "src/R/generate_asvs.r",
                "--wor_dir",
                wor_dir,
            ],
            check=True,
            # Do NOT use capture_output or stdout/stderr arguments
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

    # taxonomy_df["Taxon"] = taxonomy_df["Taxon"].str.replace("; ", ";", regex=False)
    # taxonomy_df["Taxon"] = taxonomy_df["Taxon"].replace(
    #     r";p__|;c__|;o__|;f__|;g__|;s__", ";", regex=True
    # )
    # taxonomy_df["Taxon"] = taxonomy_df["Taxon"].replace(r"d__|k__", "", regex=True)
    # taxonomy_df["Taxon"] = taxonomy_df["Taxon"].replace(
    #     to_replace="Unassigned", value="unclassified"
    # )

    return taxonomy_df


def classify_taxonomy_nb(
    rep_seq_fasta: str, num_threads: int = 1, output_path: Optional[str] = None
):
    """
    Perform taxonomic classification for 16S representative sequences.

    Args:
        rep_seq_fasta (str): Path to the FASTA file containing representative sequences.
        num_threads (int): Number of parallel jobs for classification.
        output_path (str, optional): Path to save the taxonomy classification results as CSV.

    Returns:
        pd.DataFrame: DataFrame containing taxonomy classification results.
    """
    # skip if output_path is provided and already exists
    if output_path is not None and os.path.exists(output_path):
        print(f"✅ Taxonomy classification already exists at {output_path}, skipping.")
        return pd.read_csv(output_path, index_col=0)

    # load the representative sequences as QIIME2 Artifact
    rep_seq_artifact = load_representative_sequences(rep_seq_fasta)

    # Set the classifier path for 16S
    classifier_path = "/home/maurice/resources/marker_gene_databases/16s/gtdb/226.0/v3_v4/qiime2/ssu_all_r226_extracted_derep_classifier.qza"

    # Load the classifier artifact
    print(f"Loading classifier artifact from: {classifier_path}")
    classifier_artifact = Artifact.load(classifier_path)

    # Perform taxonomic classification
    print("Classifying 16S representative sequences...")
    taxonomy_classification = feature_classifier.methods.classify_sklearn(
        classifier=classifier_artifact,
        reads=rep_seq_artifact,
        confidence=0.6,
        n_jobs=num_threads,
    )

    print("Taxonomic classification completed.")

    # Convert taxonomy to a pandas DataFrame
    taxonomy_df = taxonomy_classification.classification.view(pd.DataFrame)

    # Reformat the taxonomy DataFrame
    taxonomy_df = reformat_taxonomy(taxonomy_df)

    # Save to CSV if output_path is provided
    if output_path is not None:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        print(f"Saving taxonomy classification results to: {output_path}")
        taxonomy_df.to_csv(output_path, index=True)
        print(f"Taxonomy classification results saved to: {output_path}")

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


# def filter_fasta_sequences(wor_dir, feature_ids_to_filter):
#     """
#     Filter sequences from a FASTA file based on feature IDs to exclude.

#     Args:
#         wor_dir (str): The working directory containing the representative sequences.
#         feature_ids_to_filter (list): List of feature IDs to remove.
#     """
#     # Construct input and output FASTA paths
#     input_fasta = os.path.join(
#         wor_dir,
#         "representative_sequences",
#         "asv_sequences.fasta",
#     )
#     output_fasta = os.path.join(
#         wor_dir,
#         "representative_sequences",
#         "asv_filtered_sequences.fasta",
#     )

#     # Ensure the input FASTA file exists
#     if not Path(input_fasta).exists():
#         raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")

#     # Create the output directory if it doesn't exist
#     Path(os.path.dirname(output_fasta)).mkdir(parents=True, exist_ok=True)

#     # Open the output FASTA file and filter sequences
#     with open(output_fasta, "w") as output_handle:
#         for record in SeqIO.parse(input_fasta, "fasta"):
#             if record.id not in feature_ids_to_filter:
#                 SeqIO.write(record, output_handle, "fasta")

#     print(f"Filtered FASTA saved to {output_fasta}")
