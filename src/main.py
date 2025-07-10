import argparse
import os
import re
import subprocess
import sys
import warnings
from pathlib import Path
from typing import List

import biom
import h5py
import pandas as pd
from biom.table import Table

from python.differential_abundance_analysis import (
    create_taxa_tables,
    linda_daa_asv,
    linda_daa_picrust2,
    linda_daa_taxon,
)
from python.diversity_analysis import (
    calculate_alpha_diversity,
    generate_phylogenetic_tree_with_qiime2,
    perform_bray_curtis,
    perform_generalized_unifrac_distance,
    perform_nmds,
    perform_permanova_with_vegan,
    perform_phylo_rpca,
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


def convert_to_biom_format(asv_table_df: pd.DataFrame, output_file: str) -> Table:
    """
    Convert a DataFrame to BIOM format and save it to a file.

    Args:
        asv_table_df (pd.DataFrame): The ASV table (rows=ASVs, columns=samples).
        output_file (str): Path to save the BIOM file (.biom).

    Returns:
        biom.Table: The resulting BIOM-format table.
    """
    biom_table = Table(
        asv_table_df.values,
        observation_ids=asv_table_df.index.astype(str),
        sample_ids=asv_table_df.columns.astype(str),
    )

    # Write using HDF5 file handle
    with h5py.File(output_file, "w") as f:
        biom_table.to_hdf5(f, generated_by="convert_to_biom_format")

    return biom_table


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
    # keywords_to_filter = ["Chloroplast", "Mitochondrion"]

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

    # Read ASV table
    print(f"Reading ASV table from: {asv_table_file}")
    asv_table_df = pd.read_csv(asv_table_file, index_col=0)

    # Remove "OAC182-1", "OAC182-2","OAC182-3","OAC182-4","OAC182-5"
    asv_table_df = asv_table_df.drop(
        ["OAC182-1", "OAC182-2", "OAC182-3", "OAC182-4", "OAC182-5"], axis=1
    )

    # Remove saples that dont start with "OAC"
    asv_table_df = asv_table_df.filter(regex="^OAC", axis=1)

    # Define path to save the ASV table in BIOM format
    asv_table_biom = os.path.join(asv_table_dir, "asv_table.biom")

    biom_table = convert_to_biom_format(
        asv_table_df=asv_table_df,
        output_file=asv_table_biom,
    )

    # Save out merged data as csv
    merged_data_file = os.path.join(wor_dir, "metadata", "merged_data.csv")

    merged_df = prepare_metadata(
        asv_table_df=asv_table_df,
        wor_dir=wor_dir,
        output_path=merged_data_file,
    )

    ###############################
    # Perform PICRust2 analysis
    ################################

    # Define directory for PICRUSt2 analysis
    picrust2_dir = os.path.join(wor_dir, "picrust2_analysis")

    if os.path.exists(picrust2_dir):
        print(
            f"PICRUSt2 output directory {picrust2_dir} already exists. Skipping PICRUSt2 analysis."
        )
    else:
        command = [
            "mamba",
            "run",
            "-n",
            "picrust2",
            "picrust2_pipeline.py",
            "-s",
            rep_seq_fasta,
            "-i",
            asv_table_biom,
            "-o",
            picrust2_dir,
            "-p",
            "1",
        ]

        print(f"Running PICRUSt2 with command: {' '.join(command)}")
        # Run the command
        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode != 0:
            print("PICRUSt2 analysis failed! Exiting...")
            print(result.stderr)
            sys.exit(1)

        print("PICRUSt2 analysis completed successfully.")

    # Define pathway types and their parameters
    pathway_types = [
        {
            "name": "KO",
            "file": os.path.join(
                picrust2_dir, "KO_metagenome_out", "pred_metagenome_unstrat.tsv.gz"
            ),
        },
        {
            "name": "EC",
            "file": os.path.join(
                picrust2_dir, "EC_metagenome_out", "pred_metagenome_unstrat.tsv.gz"
            ),
        },
        {
            "name": "MetaCyc",
            "file": os.path.join(
                picrust2_dir, "pathways_out", "path_abun_unstrat.tsv.gz"
            ),
        },
    ]

    pathway_dfs = {}

    for pathway in pathway_types:
        name = pathway["name"]
        infile = pathway["file"]
        outfile = infile.replace(".tsv.gz", "_descrip.tsv.gz")

        # Add descriptions
        command = [
            "mamba",
            "run",
            "-n",
            "picrust2",
            "add_descriptions.py",
            "-i",
            infile,
            "-m",
            name,
            "-o",
            outfile,
        ]
        print(f"Adding descriptions to {name} with command: {' '.join(command)}")
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Adding descriptions to {name} failed! Exiting...")
            print(result.stderr)
            sys.exit(1)
        print(f"Descriptions added to {name} successfully.")

        # Read and re-index
        df = pd.read_csv(outfile, sep="\t", index_col=0)

        if pathway["name"] == "KO":
            df.index = [
                f"{idx.split(':')[1]} | {row['description']}"
                if ":" in idx
                else f"{idx} | {row['description']}"
                for idx, row in df.iterrows()
            ]
        elif pathway["name"] == "EC":
            df.index = df.index.astype(str) + " | " + df["description"].astype(str)
        elif pathway["name"] == "MetaCyc":
            df.index = df.index.astype(str) + " | " + df["description"].astype(str)
        else:
            raise ValueError(f"Unknown pathway type: {pathway['name']}")
        df = df.drop(columns=["description"])
        pathway_dfs[name] = df

    # Unpack for later use
    kegg_pathways_df = pathway_dfs["KO"]
    enzyme_commission_df = pathway_dfs["EC"]
    metabolic_pathways_df = pathway_dfs["MetaCyc"]

    ###############################
    # Ecological diversity analysis
    ################################

    # Define diversity metircs folder
    div_met_dir = os.path.join(wor_dir, "diversity_metrics")

    # Define directory for phlogenetics anaylsis
    phylo_dir = os.path.join(wor_dir, "phylogenetic_tree")

    rooted_tree, unrooted_tree, rooted_tree_newick_path = (
        generate_phylogenetic_tree_with_qiime2(
            input_fasta=rep_seq_fasta, threads=threads, output_dir=phylo_dir
        )
    )

    print(rooted_tree)
    print(unrooted_tree)
    print(rooted_tree_newick_path)

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

        # Define output file for the ASV table in BIOM format
        biom_file_path = os.path.join(
            asv_table_dir, f"asv_table_location_{location}.biom"
        )

        # Convert the ASV table to BIOM format
        asv_table_location_biom = convert_to_biom_format(
            asv_table_df=asv_table_location,
            output_file=biom_file_path,
        )
        #########################################
        # Differential Abundance Analysis (LINDA)
        #########################################

        # Define output file for DAA results
        daa_file = os.path.join(
            wor_dir,
            "daa_results",
            "diagnosis",
            "asv_level",
            f"daa_asv_results_location_{location}.csv",
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
                "diagnosis",
                f"{rank}_level",
                f"daa_{rank}_results_location_{location}.csv",
            )

            taxon_daa_results_df, taxon_daa_results_fil_df = linda_daa_taxon(
                taxa_table_df,
                location_df,
                primary_variable="Diagnosis",
                output_file=output_file,
            )

        # Differential Abundance Analysis (LINDA) on KEGG pathways
        print("Performing DAA on KEGG pathways...")

        kegg_table_location = kegg_pathways_df[location_df.index]

        # Define output file for DAA results
        kegg_daa_file = os.path.join(
            wor_dir,
            "daa_results",
            "diagnosis",
            "kegg_level",
            f"daa_kegg_results_location_{location}.csv",
        )

        # Run LINDA analysis on KEGG pathways
        kegg_daa_results_df, kegg_daa_results_fil_df = linda_daa_picrust2(
            function_table_df=kegg_table_location,
            metadata_df=location_df,
            primary_variable="Diagnosis",
            output_file=kegg_daa_file,
        )

        # Differential Abundance Analysis (LINDA) on EC pathways
        print("Performing DAA on EC pathways...")

        ec_table_location = enzyme_commission_df[location_df.index]

        # Define output file for DAA results
        ec_daa_file = os.path.join(
            wor_dir,
            "daa_results",
            "diagnosis",
            "ec_level",
            f"daa_ec_results_location_{location}.csv",
        )

        # Run LINDA analysis on EC pathways
        ec_daa_results_df, ec_daa_results_fil_df = linda_daa_picrust2(
            function_table_df=ec_table_location,
            metadata_df=location_df,
            primary_variable="Diagnosis",
            output_file=ec_daa_file,
        )

        # Differential Abundance Analysis (LINDA) on MetaCyc pathways
        print("Performing DAA on MetaCyc pathways...")
        meta_table_location = metabolic_pathways_df[location_df.index]

        # Define output file for DAA results
        meta_daa_file = os.path.join(
            wor_dir,
            "daa_results",
            "diagnosis",
            "metacyc_level",
            f"daa_metacyc_results_location_{location}.csv",
        )

        # Run LINDA analysis on MetaCyc pathways
        meta_daa_results_df, meta_daa_results_fil_df = linda_daa_picrust2(
            function_table_df=meta_table_location,
            metadata_df=location_df,
            primary_variable="Diagnosis",
            output_file=meta_daa_file,
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

        # Perform NMDS on the Bray-Curtis distance matrix
        coordinates_df, normalized_stress = perform_nmds(distance_matrix=bray_curtis_dm)

        # Define a title dictionary for the NMDS plot
        title_dict = {
            "primary_variable": "Diagnosis",
            "title": f"NMDS Plot of Bray-Curtis distances for location {location}",
        }

        nmds_plot = plot_nmds_with_ggplot2(
            nmds_coordinates=coordinates_df,
            metadata=location_df,
            title_dict=title_dict,
            stress_value=normalized_stress,
            permanova_results=results_dict,
        )

        beta_figures_dir = os.path.join(
            wor_dir, "figures", "beta_diversity", "bray_curtis"
        )
        # Create the directory if it doesn't exist
        os.makedirs(beta_figures_dir, exist_ok=True)

        # Define path to save the NMDS plot
        nmds_plot_file_path = os.path.join(
            beta_figures_dir, f"bray_curtis_nmds_plot_location_{location}.png"
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
            asv_table_df=asv_table_location, output_file=rpca_distance_matrix_file
        )

        # Perform PERMANOVA on the RPCA distance matrix
        rpca_results_dict = perform_permanova_with_vegan(
            metadata=location_df,
            distance_matrix=robust_aitchison_dm,
            primary_variable="Diagnosis",
        )

        # Create a title dictionary for the plot
        title_dict = {
            "primary_variable": "Diagnosis",
            "title": f"RPCA Biplot for location {location}",
        }

        rpca_plot = plot_biplot_with_ggplot2(
            ordination=ordination,
            metadata=location_df,
            title_dict=title_dict,
            permanova_results=rpca_results_dict,
            taxonomy_table=taxonomy_df,
        )

        beta_figures_dir = os.path.join(
            wor_dir, "figures", "beta_diversity", "robust_aitchison pca"
        )
        # Create the directory if it doesn't exist
        os.makedirs(beta_figures_dir, exist_ok=True)

        # Define path to save the RPCA plot
        rpca_plot_file_path = os.path.join(
            beta_figures_dir, f"rpca_biplot_location_{location}.png"
        )
        # Save the RPCA plot
        rpca_plot.save(rpca_plot_file_path, format="PNG")
        print(f"RPCA plot saved to: {rpca_plot_file_path}")

        phylo_ordination, phylo_rpca, pruned_tree, filtered_table, filtered_taxonomy = (
            perform_phylo_rpca(
                asv_table_df=asv_table_location,
                rooted_tree_newick_path=rooted_tree_newick_path,
                taxonomy_df=taxonomy_df,
                output_file=os.path.join(
                    beta_div_dir, f"phylo_rpca_location_{location}.txt"
                ),
            )
        )

        # Perform PERMANOVA on the phylo RPCA distance matrix
        phylo_rpca_results_dict = perform_permanova_with_vegan(
            metadata=location_df,
            distance_matrix=phylo_rpca,
            primary_variable="Diagnosis",
        )

        # Create a title dictionary for the phylo RPCA plot
        phylo_title_dict = {
            "primary_variable": "Diagnosis",
            "title": f"Phylo RPCA Biplot for location {location}",
        }

        phylo_rpca_plot = plot_biplot_with_ggplot2(
            ordination=phylo_ordination,
            metadata=location_df,
            title_dict=phylo_title_dict,
            permanova_results=phylo_rpca_results_dict,
            taxonomy_table=filtered_taxonomy,
        )

        beta_figures_dir = os.path.join(
            wor_dir, "figures", "beta_diversity", "phylogenetic_robust aitchison pca"
        )

        # Create the directory if it doesn't exist
        os.makedirs(beta_figures_dir, exist_ok=True)

        # Define path to save the phylo RPCA plot
        phylo_rpca_plot_file_path = os.path.join(
            beta_figures_dir, f"phylo_rpca_biplot_location_{location}.png"
        )
        # Save the phylo RPCA plot
        phylo_rpca_plot.save(phylo_rpca_plot_file_path, format="PNG")
        print(f"Phylo RPCA plot saved to: {phylo_rpca_plot_file_path}")

        gu_dm = perform_generalized_unifrac_distance(
            biom_table_path=biom_file_path,
            rooted_tree_newick_path=rooted_tree_newick_path,
        )

        # Perform PERMANOVA on the Generalized UniFrac distance matrix
        gu_results_dict = perform_permanova_with_vegan(
            metadata=location_df,
            distance_matrix=gu_dm,
            primary_variable="Diagnosis",
        )
        # Perform NMDS on the Generalized UniFrac distance matrix
        gu_coordinates_df, gu_normalized_stress = perform_nmds(distance_matrix=gu_dm)

        # Create a title dictionary for the Generalized UniFrac NMDS plot
        gu_title_dict = {
            "primary_variable": "Diagnosis",
            "title": f"NMDS Plot of Generalized UniFrac distances for location {location}",
        }

        # Plot NMDS for Generalized UniFrac
        gu_nmds_plot = plot_nmds_with_ggplot2(
            nmds_coordinates=gu_coordinates_df,
            metadata=location_df,
            title_dict=gu_title_dict,
            stress_value=gu_normalized_stress,
            permanova_results=gu_results_dict,
        )

        beta_figures_dir = os.path.join(
            wor_dir, "figures", "beta_diversity", "generalized_unifrac"
        )

        # Create the directory if it doesn't exist
        os.makedirs(beta_figures_dir, exist_ok=True)
        # Define path to save the Generalized UniFrac NMDS plot
        gu_nmds_plot_file_path = os.path.join(
            beta_figures_dir, f"generalized_unifrac_nmds_plot_location_{location}.png"
        )
        # Save the Generalized UniFrac NMDS plot
        gu_nmds_plot.save(gu_nmds_plot_file_path, format="PNG")
        print(f"Generalized UniFrac NMDS plot saved to: {gu_nmds_plot_file_path}")

    # Define directory for results
    # Copy the results to Dropboxs
    copy_to_dropbox(
        source=wor_dir,
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
