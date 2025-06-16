import argparse
import os
import re
import subprocess
import sys
import uuid
import warnings
from pathlib import Path
from typing import List

import biom
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import qiime2
import seaborn as sns
from Bio import SeqIO
from PIL import Image
from rpy2.robjects import StrVector, pandas2ri, r
from rpy2.robjects.vectors import ListVector
from skbio import DistanceMatrix, OrdinationResults

from python.differential_abundance_analysis import (
    create_taxa_tables,
    linda_daa_asv,
    linda_daa_taxon,
)
from python.diversity_analysis import (
    calculate_alpha_diversity,
    generate_phylogenetic_tree_with_qiime2,
    perform_bray_curtis,
    perform_nmds,
    perform_permanova_with_vegan,
    perform_rpca,
    rarefy_df,
)
from python.plotting import (
    plot_alpha_diversity_boxplots_with_ggplot2,
    plot_biplot_with_ggplot2,
    plot_nmds_with_ggplot2,
)
from python.processing import classify_taxonomy_nb, generate_asv_table, remove_primers
from python.qc_functions import fastq_qc


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


def prepare_metadata(
    asv_table_df: pd.DataFrame,
    wor_dir: str,
    output_path: str,
) -> pd.DataFrame:
    """
    Prepare and merge sample metadata with ASV table sample names.

    Args:
        asv_table_df (pd.DataFrame): The ASV table with sample columns.
        wor_dir (str): Working directory path.
        output_path (str, optional): Path to save the merged metadata CSV. If None, saves to default location.

    Returns:
        pd.DataFrame: Merged DataFrame with metadata and sample info.
    """

    # if output file exists, read it and return
    if os.path.exists(output_path):
        print(f"Metadata file already exists at {output_path}. Reading existing file.")
        return pd.read_csv(output_path, index_col=0)

    # Define path to metadata file
    metadata_file = os.path.join(wor_dir, "metadata", "metadata.csv")

    # Read metadata file
    metadata_df = pd.read_csv(metadata_file)
    # Rename the first column to 'patient_id'
    metadata_df.rename(columns={metadata_df.columns[0]: "patient_id"}, inplace=True)

    diagnosis_map = {
        "1": "Healthy",
        "2": "GORD",
        "3": "BO",
        "4": "Dysplasia",
        "5": "OAC",
        "6": "Metastatic",
    }

    metadata_df["Diagnosis"] = metadata_df["Diagnosis"].replace(diagnosis_map)

    # Extract sample names from the ASV table as a dataframe
    sample_data = asv_table_df.columns.to_frame(name="sample_name")
    # Get patient IDs from sample names by removing everything after the first -
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
        "BO",
        "Dysplasia",
        "OAC",
        "Metastatic",
    ]
    merged_df["Diagnosis"] = pd.Categorical(
        merged_df["Diagnosis"], categories=ordered_levels, ordered=True
    )
    # Sort the DataFrame by the Diagnosis column
    merged_df.sort_values(by="Diagnosis", inplace=True)

    merged_df.to_csv(output_path, index=True)
    print(f"Merged metadata saved to {output_path}")

    return merged_df


def depth_filtering(table: pd.DataFrame, depth_threshold: float) -> pd.DataFrame:
    """Remove samples from feature table with read depth below threshold."""
    sample_depths = table.sum(axis=0)

    # Remove samples with read depth below threshold
    return table.loc[:, sample_depths >= depth_threshold]


def prevalance_filtering(
    otu_table: pd.DataFrame, prevalence_threshold: float
) -> pd.DataFrame:
    """Remove features from feature table with prevalence below threshold."""
    feature_prevalence = (otu_table > 0).sum(axis=1) / otu_table.shape[1] * 100

    # Remove features with prevalence below threshold
    return otu_table.loc[feature_prevalence >= prevalence_threshold, :]


def convert_to_biom_format(asv_table_df: pd.DataFrame, output_file: str) -> biom.Table:
    """
    Convert a DataFrame to BIOM format and save it to a file.

    Args:
        asv_table (pd.DataFrame): The ASV table to convert.
        output_file (str): The path to save the BIOM file.
    """

    biom_table = biom.Table(
        asv_table_df.values,
        observation_ids=asv_table_df.index.astype(str),
        sample_ids=asv_table_df.columns.astype(str),
    )


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
    threads, wor_dir = parse_arguments()

    # Keywords to filter
    keywords_to_filter = ["Chloroplast", "Mitochondrion"]

    # Display the values for now (can be used later in the script)
    print(f"Number of threads: {threads}")
    print(f"Working directory: {wor_dir}")

    # Define directory  with demultiplexed FASTQ files
    demux_dir = os.path.join(wor_dir, "fastq_files", "demultiplexed")

    # Loop through each fastq file in the demux_dir and gather qc with count_nreads_fastq
    for fastq_file in Path(demux_dir).glob("*.fastq.gz"):
        base_name = fastq_file.name.replace(".fastq.gz", "")

        print(f"Processing file: {base_name}")

        # define output directory for QC results
        out_dir = os.path.join(wor_dir, "qc", "demultiplexed", base_name)

        fastq_qc(fastq_file=fastq_file, threads=threads, out_dir=out_dir)
    print("Finished quality control for demultiplexed FASTQ files.")

    # Define the directory for primer trimmed FASTQ files
    primer_trimmed_dir = os.path.join(wor_dir, "fastq_files", "primers_removed")

    # Remove primers from demultiplexed FASTQ files
    remove_primers(
        input_directory=demux_dir,
        num_threads=threads,
        output_directory=primer_trimmed_dir,
    )

    print("Finished removing primers from demultiplexed FASTQ files.")

    # Perform quality control on primer trimmed FASTQ files
    for fastq_file in Path(primer_trimmed_dir).glob("*.fastq.gz"):
        base_name = fastq_file.name.replace(".fastq.gz", "")

        print(f"Processing primer trimmed file: {base_name}")

        # Define output directory for QC results
        out_dir = os.path.join(wor_dir, "qc", "primers_removed", base_name)

        fastq_qc(fastq_file=fastq_file, threads=threads, out_dir=out_dir)
    print("Finished quality control for primer trimmed FASTQ files.")

    # # Run the R script to generate ASV table
    generate_asv_table(wor_dir=wor_dir)

    print("Finished generating ASV table.")

    # Define the path to the representative sequences fasta file
    rep_seq_fasta = os.path.join(
        wor_dir,
        "representative_sequences",
        "asv_sequences.fasta",
    )

    # Define output file path for the taxonomy classification
    output_path = os.path.join(wor_dir, "taxonomic_classification", "taxonomy_data.csv")

    # Use a pre-trained classifier
    taxonomy_df = classify_taxonomy_nb(
        rep_seq_fasta=rep_seq_fasta, output_path=output_path
    )

    # Path to ASV table
    asv_table_dir = os.path.join(wor_dir, "tables", "asv_tables")
    asv_table_file = os.path.join(asv_table_dir, "asv_table.csv")

    # Convert ASV table BIOM format
    asv_table_biom = os.path.join(asv_table_dir, "asv_table.biom")

    # Check if the ASV table BIOM file exists, if not, convert it
    # biom convert -i input.csv -o output.biom --to-hdf5 --table-type="OTU table" --header-key taxonomy

    # Read ASV table
    print(f"Reading ASV table from: {asv_table_file}")
    asv_table_df = pd.read_csv(asv_table_file, index_col=0)

    # Remove "OAC182-1", "OAC182-2","OAC182-3","OAC182-4","OAC182-5"
    asv_table_df = asv_table_df.drop(
        ["OAC182-1", "OAC182-2", "OAC182-3", "OAC182-4", "OAC182-5"], axis=1
    )

    # Remove saples that dont start with "OAC"
    asv_table_df = asv_table_df.filter(regex="^OAC", axis=1)

    # Save out merged data as csv
    merged_data_file = os.path.join(wor_dir, "metadata", "merged_data.csv")

    merged_df = prepare_metadata(
        asv_table_df=asv_table_df,
        wor_dir=wor_dir,
        output_path=merged_data_file,
    )

    ###############################
    # Ecological diversity analysis
    ################################

    # Define diversity metircs folder
    div_met_dir = os.path.join(wor_dir, "diversity_metrics")

    # Define directory for phlogenetics anaylsis
    phylo_dir = os.path.join(wor_dir, "phylogenetic_tree")

    rooted_tree, unrooted_tree = generate_phylogenetic_tree_with_qiime2(
        input_fasta=rep_seq_fasta, threads=threads, output_dir=phylo_dir
    )

    # Filter OTU table based on read depth
    asv_table_dep_fil_df = depth_filtering(asv_table_df, depth_threshold=2000)

    # rarify the ASV table
    asv_table_rarefied, rarefied_biom_table, rarefied_table_qza = rarefy_df(
        asv_table_dep_fil_df,
        with_replacement=False,
        random_seed=1088,
    )

    #######################################
    # Caculate Alpha diversity all samples
    #######################################

    # Define dir to output alpha results
    alpha_div_dir = os.path.join(div_met_dir, "alpha_diversity")

    # Define path to save alpha diversity metrics
    alpha_diversity_file = os.path.join(
        alpha_div_dir,
        "alpha_diversity_metrics.csv",
    )

    # Calculate alpha diversity metrics
    alpha_diversity_df = calculate_alpha_diversity(
        normalized_table_df=asv_table_rarefied,
        normalized_table_qza=rarefied_table_qza,
        rooted_tree=rooted_tree,
        output_file=alpha_diversity_file,
    )

    #######################################################
    # Split data by sample location
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

        #########################################
        # Differential Abundance Analysis (LINDA)
        #########################################

        # Define output file for DAA results
        daa_file = os.path.join(
            wor_dir, "daa_results", "diagnosis", f"daa_results_location_{location}.csv"
        )

        # Run LINDA analysis
        daa_results_taxa_df, daa_results_taxa_fil_df = linda_daa_asv(
            asv_table_location,
            location_df,
            primary_variable="Diagnosis",
            taxonomy_table=taxonomy_df,
            output_file=daa_file,
        )

        # Plot the DAA results
        print("Plotting DAA results...")

        # Taxonomic ranks to process
        ranks = ["species", "genus", "family", "order", "class", "phylum"]

        for rank in ranks:
            print(f"Processing taxonomic rank: {rank}")

            # Define output file path for taxa tables
            taxa_tables_file = os.path.join(
                wor_dir, "tables", "taxa_tables", f"{rank}_table.csv"
            )

            # Create taxa tables
            taxa_table_df = create_taxa_tables(
                asv_table_df,
                taxonomy_df,
                rank=rank,
                output_file=taxa_tables_file,
            )

            # Define output file for DAA results
            output_file = os.path.join(
                wor_dir,
                "daa_results",
                "taxa",
                f"daa_results_{rank}_location_{location}.csv",
            )

            taxon_daa_results_df, taxon_daa_results_fil_df = linda_daa_taxon(
                taxa_table_df,
                location_df,
                primary_variable="Diagnosis",
                output_file=output_file,
            )

        ######################
        # Alpha Diversity
        ######################

        # Subset the alpha diversity DataFrame for the current location
        alpha_diversity_location_df = alpha_diversity_df[
            alpha_diversity_df["sample_id"].isin(location_df.index)
        ]

        # Plot alpha diversity boxplots for the current location
        img = plot_alpha_diversity_boxplots_with_ggplot2(
            alpha_diversity_location_df, location_df, "Diagnosis"
        )

        # Define directory for figures
        alpha_figures_dir = os.path.join(wor_dir, "figures", "alpha_diversity")
        # Create the directory if it doesn't exist
        os.makedirs(alpha_figures_dir, exist_ok=True)

        # Save the plot as a
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
        beta_div_dir = os.path.join(wor_dir, "diversity_metrics", "beta_diversity")

        # Create the directory if it doesn't exist
        os.makedirs(beta_div_dir, exist_ok=True)

        # rarify the ASV table
        rarefied_df, rarefied_biom_table, rarefied_table_qza = rarefy_df(
            asv_table_location,
            with_replacement=False,
            random_seed=1088,
        )

        output_file = os.path.join(
            beta_div_dir, f"bray_curtis_distance_matrix_{location}.txt"
        )

        # Perform Bray-Curtis distance calculation
        bray_curtis_dm = perform_bray_curtis(
            rarefied_table_df=rarefied_df, output_file=output_file
        )

        results_dict = perform_permanova_with_vegan(
            metadata=location_df,
            distance_matrix=bray_curtis_dm,
            primary_variable="Diagnosis",
        )

        coordinates_df, normalized_stress = perform_nmds(distance_matrix=bray_curtis_dm)

        nmds_plot = plot_nmds_with_ggplot2(
            nmds_coordinates=coordinates_df,
            metadata=location_df,
            primary_variable="Diagnosis",
            stress_value=normalized_stress,
            permanova_results=results_dict,
        )

        beta_figures_dir = os.path.join(wor_dir, "figures", "beta_diversity")
        # Create the directory if it doesn't exist
        os.makedirs(beta_figures_dir, exist_ok=True)

        # Define path to save the NMDS plot
        nmds_plot_file_path = os.path.join(
            beta_figures_dir, f"nmds_plot_location_{location}.png"
        )
        # Save the NMDS plot
        nmds_plot.save(nmds_plot_file_path, format="PNG")
        print(f"NMDS plot saved to: {nmds_plot_file_path}")

        # define outfile for RPCA distance matrix
        rpca_distance_matrix_file = os.path.join(
            beta_div_dir, f"rpca_location_{location}_distance_matrix.txt"
        )

        # Perform RPCA
        ordination, robust_aitchison_dm = perform_rpca(
            asv_table_df=rarefied_df, output_file=rpca_distance_matrix_file
        )

        # Perform PERMANOVA on the RPCA distance matrix
        rpca_results_dict = perform_permanova_with_vegan(
            metadata=location_df,
            distance_matrix=robust_aitchison_dm,
            primary_variable="Diagnosis",
        )

        rpca_plot = plot_biplot_with_ggplot2(
            ordination=ordination,
            metadata=location_df,
            primary_variable="Diagnosis",
            permanova_results=rpca_results_dict,
            taxonomy_table=taxonomy_df,
        )

        # Define path to save the RPCA plot
        rpca_plot_file_path = os.path.join(
            beta_figures_dir, f"rpca_plot_location_{location}.png"
        )
        # Save the RPCA plot
        rpca_plot.save(rpca_plot_file_path, format="PNG")
        print(f"RPCA plot saved to: {rpca_plot_file_path}")

    # Copy the results to Dropboxs
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
