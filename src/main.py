import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import List

import h5py  # type: ignore
import pandas as pd # type: ignore
from biom.table import Table  # type: ignore

from python.beta_diversity_analysis import (
    calculate_bray_curtis_distance,
    calculate_generalized_unifrac_distance,
    calculate_jaccard_distance,
    calculate_unweighted_unifrac_distance,
    perform_pairwise_comparisons_permanova,
    perform_permanova_with_vegan,
    perform_phylogenetic_rpca,
    perform_rpca,
)
from python.differential_abundance_analysis import (
    create_taxa_tables,
    linda_daa_asv,
    linda_daa_picrust2,
    linda_daa_taxon,
)
from python.diversity_analysis import (
    calculate_alpha_diversity,
    convert_dataframe_to_qiime2_artifact,
    generate_phylogenetic_tree_with_qiime2,
    paired_alpha_diversity_calculations,
    perform_ctf,
    perform_nmds,
    rarefy_df,
    unpaired_alpha_diversity_calculations,
)
from python.plotting import (
    plot_biplot_with_ggplot2,
    plot_heatmap_of_alpha_diversity,
    plot_heatmap_of_daa,
    plot_nmds_with_ggplot2,
)
from python.processing import classify_taxonomy_nb, generate_asv_table, remove_primers
from python.qc_functions import fastq_qc, plot_read_counts_primer_trimmed


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
    sample_data = asv_table_df.columns.to_frame(
        name="sample_id"
    )  # Changed from sample_name
    # Get patient IDs from sample names by removing everything after the first -
    sample_data["patient_id"] = sample_data["sample_id"].str.extract(
        r"^([^-]+)"
    )  # Changed from sample_name
    # Get sample locations from sample names by removing everything before the first -
    sample_data["sample_location"] = sample_data["sample_id"].str.extract(
        r"-(.*)"
    )  # Changed from sample_name

    # Merge sample data with metadata
    merged_df = pd.merge(
        sample_data,
        metadata_df,
        left_on="patient_id",
        right_on="patient_id",
        how="inner",
    )

    # Set sample_id as the index
    merged_df.set_index("sample_id", inplace=True)  # Changed from sample_name

    # In location column, append "biopsy_location_" before each location number
    merged_df["sample_location"] = "biopsy_location_" + merged_df[
        "sample_location"
    ].astype(str)

    # Set order level of locations
    location_order = [
        "biopsy_location_1",
        "biopsy_location_2",
        "biopsy_location_3",
        "biopsy_location_4",
        "biopsy_location_5",
    ]

    merged_df["sample_location"] = pd.Categorical(
        merged_df["sample_location"], categories=location_order, ordered=True
    )

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
    # Order the dataframe by location and then diagnosis
    merged_df.sort_values(by=["sample_location", "Diagnosis"], inplace=True)

    # Convert individual_id_column and state_column to string type
    merged_df["Diagnosis"] = merged_df["Diagnosis"].astype(str)
    merged_df["sample_location"] = merged_df["sample_location"].astype(str)

    merged_df.to_csv(output_path, index=True)
    print(f"Merged metadata saved to {output_path}")

    return merged_df


def depth_filtering(table: pd.DataFrame, depth_threshold: float) -> pd.DataFrame:
    """Remove samples from feature table with read depth below threshold."""
    sample_depths = table.sum(axis=0)

    # Remove samples with read depth below threshold
    return table.loc[:, sample_depths >= depth_threshold]


def convert_to_biom_format(asv_table_df: pd.DataFrame, output_file: str) -> None:
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

    ###############################################################################
    ###############################################################################
    ###############################################################################
    # 1 Raw data processing and data generation
    ################################################################################
    ################################################################################
    ################################################################################

    #################################
    #################################
    # 1.1 Primer removal and QC
    #################################
    #################################

    ####################################################
    # 1.1.1 Quality control on demultiplexed FASTQ files
    ####################################################

    # Check if demultiplexed QC directory exists and has files
    demux_qc_dir = os.path.join(wor_dir, "qc", "demultiplexed")

    if os.path.exists(demux_qc_dir):
        existing_qc_files = list(Path(demux_qc_dir).glob("*"))
        if existing_qc_files:
            print(
                f"Demultiplexed QC directory already contains {len(existing_qc_files)} files. Skipping demultiplexed QC."
            )
        else:
            print(
                "Demultiplexed QC directory exists but is empty. Running demultiplexed QC..."
            )
            # Loop through each fastq file in the demux_dir and gather qc with count_nreads_fastq
            for fastq_file in Path(demux_dir).glob("*.fastq.gz"):
                base_name = fastq_file.name.replace(".fastq.gz", "")

                print(f"Processing file: {base_name}")

                # define output directory for QC results
                out_dir = os.path.join(wor_dir, "qc", "demultiplexed", base_name)

                fastq_qc(fastq_file=fastq_file, threads=threads, out_dir=out_dir)
            print("Finished quality control for demultiplexed FASTQ files.")
    else:
        print("Demultiplexed QC directory doesn't exist. Running demultiplexed QC...")
        # Loop through each fastq file in the demux_dir and gather qc with count_nreads_fastq
        for fastq_file in Path(demux_dir).glob("*.fastq.gz"):
            base_name = fastq_file.name.replace(".fastq.gz", "")

            print(f"Processing file: {base_name}")

            # define output directory for QC results
            out_dir = os.path.join(wor_dir, "qc", "demultiplexed", base_name)

            fastq_qc(fastq_file=fastq_file, threads=threads, out_dir=out_dir)
        print("Finished quality control for demultiplexed FASTQ files.")

    #####################################################
    # 1.1.2 Primer removal from demultiplexed FASTQ files
    ######################################################

    # Define the directory for primer trimmed FASTQ files
    primer_trimmed_dir = os.path.join(wor_dir, "fastq_files", "primers_removed")

    # Check if primer trimmed directory exists and has files
    if os.path.exists(primer_trimmed_dir):
        existing_files = list(Path(primer_trimmed_dir).glob("*.fastq.gz"))
        if existing_files:
            print(
                f"Primer trimmed directory already contains {len(existing_files)} files. Skipping primer removal."
            )
        else:
            print(
                "Primer trimmed directory exists but is empty. Running primer removal..."
            )
            remove_primers(
                input_directory=demux_dir,
                num_threads=threads,
                output_directory=primer_trimmed_dir,
            )
            print("Finished removing primers from demultiplexed FASTQ files.")
    else:
        print("Primer trimmed directory doesn't exist. Running primer removal...")
        remove_primers(
            input_directory=demux_dir,
            num_threads=threads,
            output_directory=primer_trimmed_dir,
        )
        print("Finished removing primers from demultiplexed FASTQ files.")

    #####################################################
    # 1.1.3 Quality control on primer trimmed FASTQ files
    #####################################################

    # Perform quality control on primer trimmed FASTQ files
    primer_trimmed_qc_dir = os.path.join(wor_dir, "qc", "primers_removed")

    # Check if primer trimmed QC directory exists and has files
    if os.path.exists(primer_trimmed_qc_dir):
        existing_qc_files = list(Path(primer_trimmed_qc_dir).glob("*"))
        if existing_qc_files:
            print(
                f"Primer trimmed QC directory already contains {len(existing_qc_files)} files. Skipping primer trimmed QC."
            )
        else:
            print(
                "Primer trimmed QC directory exists but is empty. Running primer trimmed QC..."
            )
            for fastq_file in Path(primer_trimmed_dir).glob("*.fastq.gz"):
                base_name = fastq_file.name.replace(".fastq.gz", "")
                print(f"Processing primer trimmed file: {base_name}")
                out_dir = os.path.join(wor_dir, "qc", "primers_removed", base_name)
                fastq_qc(fastq_file=fastq_file, threads=threads, out_dir=out_dir)
            print("Finished quality control for primer trimmed FASTQ files.")
    else:
        print("Primer trimmed QC directory doesn't exist. Running primer trimmed QC...")
        for fastq_file in Path(primer_trimmed_dir).glob("*.fastq.gz"):
            base_name = fastq_file.name.replace(".fastq.gz", "")
            print(f"Processing primer trimmed file: {base_name}")
            out_dir = os.path.join(wor_dir, "qc", "primers_removed", base_name)
            fastq_qc(fastq_file=fastq_file, threads=threads, out_dir=out_dir)
        print("Finished quality control for primer trimmed FASTQ files.")

    #########################################################
    # 1.1.4 Plot read counts before and after primer trimming
    ##########################################################

    # Plot read counts before and after primer trimming
    plot_read_counts_primer_trimmed(
        folder1=os.path.join(wor_dir, "qc", "demultiplexed"),
        folder2=os.path.join(wor_dir, "qc", "primers_removed"),
        output_file=os.path.join(wor_dir, "qc", "read_counts_primer_trimmed.csv"),
        plot_file=os.path.join(wor_dir, "qc", "read_counts_primer_trimmed.png"),
    )

    #########################################################
    #########################################################
    # 1.2 ASV table generation with DADA2
    #########################################################
    #########################################################

    # # Run the R script to generate ASV table
    generate_asv_table(wor_dir=wor_dir)

    print("Finished generating ASV table.")

    #######################################################
    #######################################################
    # 1.3 Prepare ASV table and merge with metadata
    #######################################################
    #######################################################

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

    asv_table_dep_fil_df = depth_filtering(asv_table_df, depth_threshold=2000)

    # Define path to save the ASV table in BIOM format
    asv_table_biom = os.path.join(asv_table_dir, "asv_table.biom")

    convert_to_biom_format(
        asv_table_df=asv_table_dep_fil_df,
        output_file=asv_table_biom,
    )

    # Save out merged data as csv
    merged_data_file = os.path.join(wor_dir, "metadata", "merged_data.csv")

    merged_df = prepare_metadata(
        asv_table_df=asv_table_dep_fil_df,
        wor_dir=wor_dir,
        output_path=merged_data_file,
    )

    # Subset asv_table_dep_fil_df to only include samples in merged_df
    asv_table_dep_fil_df = asv_table_dep_fil_df[merged_df.index]

    ##########################################################################
    ##########################################################################
    ##########################################################################
    # 2. Taxonomy classification, functional prediction and phylogenetics
    ###########################################################################
    ###########################################################################

    #####################################################
    #####################################################
    # 2.1 Taxonomy classification
    #####################################################
    #####################################################

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

    #######################################################
    #######################################################
    # 2.5 PICRust2 analysis
    #######################################################
    #######################################################

    #######################################
    # 2.5.1 Run PICRUSt2 analysis
    #######################################

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

        print("Running PICRUSt2...")
        # Run the command
        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode != 0:
            print("PICRUSt2 analysis failed! Exiting...")
            print(result.stderr)
            sys.exit(1)

        print("PICRUSt2 analysis completed successfully.")

    ###########################################
    # 2.5.2 Add descriptions to PICRUSt2 output
    ###########################################

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

        # Skip if output file already exists
        if os.path.exists(outfile):
            print(
                f"Description file already exists. Skipping description addition for {name}."
            )
        else:
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

        ###########################################
        # 2.5.3 Read and re-index PICRUSt2 output
        ###########################################

        # Read and re-index (always do this regardless of whether we skipped or ran the command)
        df = pd.read_csv(outfile, sep="\t", index_col=0)

        if pathway["name"] == "KO":
            df.index = [
                (
                    f"{idx.split(':')[1]} | {row['description']}"
                    if ":" in idx
                    else f"{idx} | {row['description']}"
                )
                for idx, row in df.iterrows()
            ]
        elif pathway["name"] == "EC":
            df.index = df.index.astype(str) + " | " + df["description"].astype(str)
        elif pathway["name"] == "MetaCyc":
            df.index = df["description"].astype(str)
        else:
            raise ValueError(f"Unknown pathway type: {pathway['name']}")
        df = df.drop(columns=["description"])
        pathway_dfs[name] = df

    # Extract individual DataFrames
    kegg_pathways_df = pathway_dfs["KO"]
    enzyme_commission_df = pathway_dfs["EC"]
    metabolic_pathways_df = pathway_dfs["MetaCyc"]

    ###############################
    # 2.6 Phylogenetic tree
    ################################

    # Define diversity metrics folder
    div_met_dir = os.path.join(wor_dir, "diversity_metrics")

    # Define directory for phylogenetics analysis
    phylo_dir = os.path.join(wor_dir, "phylogenetic_tree")

    rooted_tree, unrooted_tree, rooted_tree_newick_path = (
        generate_phylogenetic_tree_with_qiime2(
            input_fasta=rep_seq_fasta, threads=threads, output_dir=phylo_dir
        )
    )

    # rarify the ASV table
    asv_table_rarefied = rarefy_df(
        asv_table_dep_fil_df,
        with_replacement=False,
        random_seed=1088,
    )

    # Convert the rarefied ASV table to a QIIME2 artifact
    rarefied_table_qza = convert_dataframe_to_qiime2_artifact(asv_table_rarefied)

    #######################################
    # 1.7 Alpha diversity metrics
    #######################################

    # Define dir to output alpha results
    alpha_div_dir = os.path.join(div_met_dir, "alpha_diversity")

    # Define path to save alpha diversity metrics
    alpha_diversity_file = os.path.join(
        alpha_div_dir,
        "alpha_diversity_metrics.csv",
    )

    if os.path.exists(alpha_diversity_file):
        print(
            f"Alpha diversity metrics file already exists at {alpha_diversity_file}. Skipping alpha diversity calculations."
        )
        alpha_diversity_df = pd.read_csv(alpha_diversity_file)
    else:
        print("Calculating alpha diversity metrics...")

        # Calculate alpha diversity metrics
        alpha_diversity_df = calculate_alpha_diversity(
            normalized_table_df=asv_table_rarefied,
            normalized_table_qza=rarefied_table_qza,
            rooted_tree=rooted_tree,
            output_file=alpha_diversity_file,
        )
        print(f"Alpha diversity metrics saved to {alpha_diversity_file}")

    #######################################################################
    #######################################################################
    #######################################################################
    # 3. Data analysis
    #######################################################################
    #######################################################################

    #######################################################################
    #######################################################################
    # 3.1 All samples together analysis
    #######################################################################
    #######################################################################

    ######################################################################
    # 3.1.1 Beta diversity analysis of all samples together
    #######################################################

    # Define beta diversity directory
    beta_div_dir = os.path.join(div_met_dir, "beta_diversity")

    #########################
    # 3.1.1.1 CTF analysis
    #########################
    # Define dir to output directory for CTF results
    ctf_output_dir = os.path.join(beta_div_dir, "ctf")

    # Path to figures directory
    fig_dir = os.path.join(wor_dir, "figures")

    # Path to beta diversity figures directory
    beta_fig_dir = os.path.join(fig_dir, "beta_diversity")

    # Path to CTF biplot directory
    ctf_biplot_dir = os.path.join(beta_fig_dir, "ctf")
    ctf_biplot_output_file = os.path.join(ctf_biplot_dir, "ctf_biplot.pdf")

    # Check if CTF biplot already exists
    if os.path.exists(ctf_biplot_output_file):
        print(
            f"CTF biplot already exists at {ctf_biplot_output_file}. Skipping CTF analysis."
        )
    else:
        # Perform CTF
        ordination_results, ctf_distance_df = perform_ctf(
            asv_count_table_df=asv_table_dep_fil_df,
            metadata_df=merged_df,
            individual_id_column="patient_id",
            state_column="sample_location",
            output_dir=ctf_output_dir,
        )

        # Perform pairwise comparisons with vegan on the CTF distance matrix
        perform_pairwise_comparisons_permanova(
            metadata=merged_df,
            distance_df=ctf_distance_df,
            primary_variable="Diagnosis",
            output_file=os.path.join(ctf_output_dir, "pairwise_permanova_results.csv"),
            all=True,
        )

        # Perform permanova  with vegan on the CTF distance matrix
        permanova_results = perform_permanova_with_vegan(
            metadata=merged_df,
            distance_df=ctf_distance_df,
            primary_variable="Diagnosis",
            strata="True",
        )

        # Plot CTF biplot
        plot_biplot_with_ggplot2(
            ordination_results=ordination_results,
            metadata=merged_df,
            title_dict={
                "primary_variable": "Diagnosis",
                "title": "CTF Biplot of Locations",
            },
            permanova_results=permanova_results,
            taxonomy_table=taxonomy_df,
            output_file=ctf_biplot_output_file,
            all=True,
        )
        print(f"CTF biplot saved to {ctf_biplot_output_file}")

    ##########################################################
    ##########################################################
    # 3.2 Comparative analysis between locations per diagnosis
    ##########################################################
    ##########################################################

    diagnoses = merged_df["Diagnosis"].unique()

    for diagnosis in diagnoses:
        ########################################################################
        # 3.2.1 Generate diagnosis specific ASV table in pandas and BIOM format
        ########################################################################
        print(f"Processing diagnosis: {diagnosis}")

        # Filter the merged DataFrame for the current diagnosis
        diagnosis_df = merged_df[merged_df["Diagnosis"] == diagnosis]

        # Filter the ASV table for the current diagnosis
        asv_table_diagnosis = asv_table_dep_fil_df[diagnosis_df.index]

        # rarify the ASV table
        rarefied_df = rarefy_df(
            asv_table_diagnosis,
            with_replacement=False,
            random_seed=1088,
        )

        # Define output file for the ASV table in BIOM format
        biom_file_path = os.path.join(
            asv_table_dir, f"asv_rarefied_table_{diagnosis}.biom"
        )

        # Convert to BIOM format
        convert_to_biom_format(
            asv_table_df=rarefied_df,
            output_file=biom_file_path,
        )
        print(f"BIOM file saved to {biom_file_path}")

        ############################################################
        # 3.2.2 Alpha diversity between locations per diagnosis
        ############################################################

        # Path to alpha divesity figures directory
        alpha_fig_dir = os.path.join(fig_dir, "alpha_diversity")

        # Define to diagnosis specific alpha diversity directory
        diagnosis_alpha_fig_dir = os.path.join(alpha_fig_dir, "location")

        # Path to file to save alpha diversity figures
        alpha_diversity_figure_file = os.path.join(
            diagnosis_alpha_fig_dir,
            f"alpha_diversity_metrics_location_comparison_{diagnosis}.pdf",
        )

        # Subset the alpha diversity DataFrame for the current diagnosis
        alpha_diversity_diagnosis_df = alpha_diversity_df[
            alpha_diversity_df["sample_id"].isin(diagnosis_df.index)
        ]

        # Path results of statiscal comparisons of alpha diversity between locations per diagnosis
        alpha_diversity_diagnosis_file = os.path.join(
            alpha_div_dir,
            "results",
            "location",
            f"alpha_diversity_location_comparison_{diagnosis}.csv",
        )

        # Calculate differences in alpha diversity between locations
        alpha_diversity_results_df = paired_alpha_diversity_calculations(
            alpha_diversity_df=alpha_diversity_diagnosis_df,
            metadata_df=diagnosis_df,
            primary_variable="sample_location",
            output_file=alpha_diversity_diagnosis_file,
        )

        plot_heatmap_of_alpha_diversity(
            data_df=alpha_diversity_results_df,
            title_dict={
                "primary_variable": "sample_location",
                "title": f"Alpha Diversity Comparison between Biopsy Locations for {diagnosis} Patients",
            },
            output_file=alpha_diversity_figure_file,
        )

        # Save the figure
        print(f"Alpha diversity figure saved to {alpha_diversity_figure_file}")

        ######################################################
        # 3.2.3 Beta diversity between locations per diagnosis
        #######################################################

        # Path to directory to output beta diversity results locations per diagnosis
        beta_div_location_dir = os.path.join(beta_div_dir, "location")

        # Path to directory to figures of beta diversity results locations per diagnosis
        beta_fig_location_dir = os.path.join(beta_fig_dir, "location")

        ################################
        # 3.2.3.1 Jaccard
        ################################

        # Path to Jaccard directory
        jaccard_dir = os.path.join(beta_div_location_dir, "jaccard")

        # Path to figures directory for Jaccard results comparing locations per diagnosis
        jaccard_fig_location_dir = os.path.join(beta_fig_location_dir, "jaccard")

        # Path to NMDS plot file comparing locations per diagnosis
        nmds_plot_file_path = os.path.join(
            jaccard_fig_location_dir,
            f"jaccard_nmds_plot_comparing_locations_{diagnosis}.pdf",
        )

        if os.path.exists(nmds_plot_file_path):
            print(
                f"Jaccard NMDS plot already exists at {nmds_plot_file_path}. Skipping Jaccard analysis for {diagnosis}."
            )
        else:
            print(f"Calculating Jaccard distances foqr {diagnosis}...")
            # Path to Jaccard distance matrix comparing locations per diagnosis
            jaccard_distance_matrix_file = os.path.join(
                jaccard_dir,
                f"jaccard_distance_matrix_comparing_locations_{diagnosis}.csv",
            )

            # Calculate Jaccard distance
            jaccard_df, jaccard_dm = calculate_jaccard_distance(
                rarefied_table_df=rarefied_df,
                output_file=jaccard_distance_matrix_file,
            )

            # Perform pairwise comparisons with vegan on the Jaccard distance matrix
            perform_pairwise_comparisons_permanova(
                metadata=diagnosis_df,
                distance_df=jaccard_df,
                primary_variable="sample_location",
                output_file=os.path.join(
                    jaccard_dir,
                    f"jaccard_pairwise_comparisons_between_locations_{diagnosis}.csv",
                ),
            )

            # Perform permanova  with vegan on the Jaccard distance matrix
            permanova_results_dict = perform_permanova_with_vegan(
                metadata=diagnosis_df,
                distance_df=jaccard_df,
                primary_variable="sample_location",
            )

            # Perform NMDS on the Jaccard distance matrix
            coordinates_df, normalized_stress = perform_nmds(distance_matrix=jaccard_dm)

            # Define a title dictionary for the NMDS plot
            title_dict = {
                "primary_variable": "sample_location",
                "title": f"NMDS Plot of Jaccard distances comparing locations in {diagnosis} Patients",
                "stress": normalized_stress,
            }

            plot_nmds_with_ggplot2(
                nmds_coordinates=coordinates_df,
                metadata=diagnosis_df,
                permanova_results=permanova_results_dict,
                title_dict=title_dict,
                output_file=nmds_plot_file_path,
            )

        ########################
        # 3.2.3.2 Bray-Curtis
        ########################

        # Path to Bray-Curtis directory
        bray_curtis_location_dir = os.path.join(beta_div_location_dir, "bray_curtis")

        # Path to figures directory for Bray-Curtis results comparing locations per diagnosis
        bray_curtis_fig_location_dir = os.path.join(
            beta_fig_location_dir, "bray_curtis"
        )

        # Path to NMDS plot file comparing locations per diagnosis
        nmds_plot_file_path = os.path.join(
            bray_curtis_fig_location_dir,
            f"bray_curtis_nmds_plot_comparing_locations_{diagnosis}.pdf",
        )

        if os.path.exists(nmds_plot_file_path):
            print(
                f"Bray-Curtis NMDS plot already exists at {nmds_plot_file_path}. Skipping Bray-Curtis analysis for {diagnosis}."
            )
        else:
            print(f"Calculating Bray-Curtis distances for {diagnosis}...")
            # Path to Bray-Curtis distance matrix comparing locations per diagnosis
            bray_curtis_distance_matrix_file = os.path.join(
                bray_curtis_location_dir,
                f"bray_curtis_distance_matrix_comparing_locations_{diagnosis}.csv",
            )

            # Perform Bray-Curtis distance calculation
            bray_curtis_df, bray_curtis_dm = calculate_bray_curtis_distance(
                rarefied_table_df=rarefied_df,
                output_file=bray_curtis_distance_matrix_file,
            )

            results_dict = perform_permanova_with_vegan(
                metadata=diagnosis_df,
                distance_df=bray_curtis_df,
                primary_variable="sample_location",
            )

            # Perform pairwise comparisons with vegan on the Bray-Curtis distance matrix
            perform_pairwise_comparisons_permanova(
                metadata=diagnosis_df,
                distance_df=bray_curtis_df,
                primary_variable="sample_location",
                output_file=os.path.join(
                    bray_curtis_location_dir,
                    f"bray_curtis_pairwise_comparisons_comparing_locations_{diagnosis}.csv",
                ),
            )

            # Perform NMDS on the Bray-Curtis distance matrix
            coordinates_df, normalized_stress = perform_nmds(
                distance_matrix=bray_curtis_dm
            )

            # Define a title dictionary for the NMDS plot
            title_dict = {
                "primary_variable": "sample_location",
                "title": f"NMDS Plot of Bray-Curtis distances comparing locations in {diagnosis} Patients",
                "stress": normalized_stress,
            }

            plot_nmds_with_ggplot2(
                nmds_coordinates=coordinates_df,
                metadata=diagnosis_df,
                title_dict=title_dict,
                permanova_results=results_dict,
                output_file=nmds_plot_file_path,
            )

        ###############################
        # 3.2.3.3 Robust Aitchison PCA
        ###############################

        # Path to Bray-Curtis directory
        rpca_location_dir = os.path.join(beta_div_location_dir, "rpca")

        # Path to figures directory for RPCA results comparing locations per diagnosis
        rpca_fig_location_dir = os.path.join(beta_fig_location_dir, "rpca")

        # Path to RPCA biplot file comparing locations per diagnosis
        rpca_biplot_file = os.path.join(
            rpca_fig_location_dir, f"rpca_biplot_comparing_locations_{diagnosis}.pdf"
        )

        if os.path.exists(rpca_biplot_file):
            print(
                f"RPCA biplot already exists at {rpca_biplot_file}. Skipping RPCA analysis for {diagnosis}."
            )
        else:
            # Path to RPCA distance matrix comparing locations per diagnosis
            rpca_distance_matrix_file = os.path.join(
                rpca_location_dir,
                f"rpca_distance_matrix_comparing_locations_{diagnosis}.csv",
            )

            # Perform RPCA
            rpca_ordination_results, robust_aitchison_distance_df = perform_rpca(
                asv_table_df=asv_table_diagnosis,
                output_file=rpca_distance_matrix_file,
            )

            # Perform pairwise comparisons with vegan on the RPCA distance matrix
            perform_pairwise_comparisons_permanova(
                metadata=diagnosis_df,
                distance_df=robust_aitchison_distance_df,
                primary_variable="sample_location",
                output_file=os.path.join(
                    rpca_location_dir,
                    f"rpca_pairwise_comparisons_between_locations_{diagnosis}.csv",
                ),
            )

            # Perform permanova  with vegan on the RPCA distance matrix
            permanova_results = perform_permanova_with_vegan(
                metadata=diagnosis_df,
                distance_df=robust_aitchison_distance_df,
                primary_variable="sample_location",
            )

            # Create dict for title and primary variable
            title_dict = {
                "primary_variable": "sample_location",
                "title": f"RPCA Biplot comparing locations for {diagnosis} Patients",
            }

            # Plot RPCA biplot
            plot_biplot_with_ggplot2(
                ordination_results=rpca_ordination_results,
                metadata=diagnosis_df,
                permanova_results=permanova_results,
                taxonomy_table=taxonomy_df,
                title_dict=title_dict,
                output_file=rpca_biplot_file,
            )

        ############################################
        # 3.2.3.4  Phylogenetic robust Aitchison PCA
        ############################################

        # Path to Phylogenetic RPCA directory
        phylo_rpca_location_dir = os.path.join(beta_div_location_dir, "phylo_rpca")

        # Path to figures directory for Phylogenetic RPCA results comparing locations per diagnosis
        phylo_rpca_fig_location_dir = os.path.join(beta_fig_location_dir, "phylo_rpca")

        # Path to Phylogenetic RPCA biplot file comparing locations per diagnosis
        phylo_rpca_biplot_file = os.path.join(
            phylo_rpca_fig_location_dir,
            f"phylo_rpca_biplot_comparing_locations_{diagnosis}.pdf",
        )

        if os.path.exists(phylo_rpca_biplot_file):
            print(
                f"Phylogenetic RPCA biplot already exists at {phylo_rpca_biplot_file}. Skipping Phylogenetic RPCA analysis for {diagnosis}."
            )
        else:
            print(f"Calculating Phylogenetic RPCA for {diagnosis}...")

            # Path to Phylogenetic RPCA distance matrix comparing locations per diagnosis
            phylo_rpca_distance_matrix_file = os.path.join(
                phylo_rpca_location_dir,
                f"phylo_rpca_distance_matrix_comparing_locations_{diagnosis}.csv",
            )

            # Perform Phylogenetic RPCA

            (
                phylo_rpca_ordination,
                phylo_rpca_distance_df,
                pruned_tree,
                filtered_table,
                filtered_taxonomy,
            ) = perform_phylogenetic_rpca(
                asv_table_df=asv_table_diagnosis,
                taxonomy_df=taxonomy_df,
                rooted_tree_newick_path=rooted_tree_newick_path,
                output_file=phylo_rpca_distance_matrix_file,
            )

            # Perform pairwise comparisons with vegan on the Phylogenetic RPCA distance matrix
            perform_pairwise_comparisons_permanova(
                metadata=diagnosis_df,
                distance_df=phylo_rpca_distance_df,
                primary_variable="sample_location",
                output_file=os.path.join(
                    phylo_rpca_location_dir,
                    f"phylo_rpca_pairwise_comparisons_between_locations_{diagnosis}.csv",
                ),
            )

            # Perform permanova  with vegan on the Phylogenetic RPCA distance matrix
            permanova_results = perform_permanova_with_vegan(
                metadata=diagnosis_df,
                distance_df=phylo_rpca_distance_df,
                primary_variable="sample_location",
            )

            # Create dict for title and primary variable
            title_dict = {
                "primary_variable": "sample_location",
                "title": f"Phylogenetic RPCA Biplot comparing locations for {diagnosis} Patients",
            }

            # Plot Phylogenetic RPCA biplot
            plot_biplot_with_ggplot2(
                ordination_results=phylo_rpca_ordination,
                metadata=diagnosis_df,
                permanova_results=permanova_results,
                taxonomy_table=taxonomy_df,
                title_dict=title_dict,
                output_file=phylo_rpca_biplot_file,
            )

        ##############################
        # 3.2.3.5  Unweighted UniFrac
        ##############################

        # Path to Unweighted UniFrac directory
        unweighted_unifrac_location_dir = os.path.join(
            beta_div_location_dir, "unweighted_unifrac"
        )

        # Path to figures directory for Unweighted UniFrac results comparing locations per diagnosis
        unweighted_unifrac_fig_location_dir = os.path.join(
            beta_fig_location_dir, "unweighted_unifrac"
        )

        if os.path.exists(nmds_plot_file_path):
            print(
                f"Unweighted UniFrac NMDS plot already exists at {nmds_plot_file_path}. Skipping Unweighted UniFrac analysis for {diagnosis}."
            )
        else:
            print(f"Calculating Unweighted UniFrac distances for {diagnosis}...")

            # Path to NMDS plot file comparing locations per diagnosis
            nmds_plot_file_path = os.path.join(
                unweighted_unifrac_fig_location_dir,
                f"unweighted_unifrac_nmds_plot_comparing_locations_{diagnosis}.pdf",
            )

            # Path to Unweighted UniFrac distance matrix comparing locations per diagnosis
            unweighted_unifrac_distance_matrix_file = os.path.join(
                unweighted_unifrac_location_dir,
                f"unweighted_unifrac_distance_matrix_comparing_locations_{diagnosis}.csv",
            )

            # Calculate Unweighted UniFrac distance
            unweighted_unifrac_df, unweighted_unifrac_dm = (
                calculate_unweighted_unifrac_distance(
                    biom_table_path=biom_file_path,
                    rooted_tree_newick_path=rooted_tree_newick_path,
                    output_file=unweighted_unifrac_distance_matrix_file,
                )
            )

            # Perform NMDS
            coordinates_df, normalized_stress = perform_nmds(
                distance_matrix=unweighted_unifrac_dm
            )

            # Perform pairwise comparisons with vegan on the Unweighted UniFrac distance matrix
            perform_pairwise_comparisons_permanova(
                metadata=diagnosis_df,
                distance_df=unweighted_unifrac_df,
                primary_variable="sample_location",
                output_file=os.path.join(
                    unweighted_unifrac_location_dir,
                    f"unweighted_unifrac_pairwise_comparisons_between_locations_{diagnosis}.csv",
                ),
            )

            # Perform permanova  with vegan on the Unweighted UniFrac distance matrix
            results_dict = perform_permanova_with_vegan(
                metadata=diagnosis_df,
                distance_df=unweighted_unifrac_df,
                primary_variable="sample_location",
            )

            # Create NMDS plot
            title_dict = {
                "primary_variable": "sample_location",
                "title": f"NMDS Plot of Unweighted UniFrac distances comparing locations in {diagnosis} Patients",
                "stress": normalized_stress,
            }

            # Plot NMDS
            plot_nmds_with_ggplot2(
                nmds_coordinates=coordinates_df,
                metadata=diagnosis_df,
                permanova_results=results_dict,
                title_dict=title_dict,
                output_file=nmds_plot_file_path,
            )

        ###############################
        # 3.2.3.5  Generalized UniFrac
        ##############################

        # Path to Generalized UniFrac directory
        generalized_unifrac_location_dir = os.path.join(
            beta_div_location_dir, "generalized_unifrac"
        )

        # Path to figures directory for Generalized UniFrac results comparing locations per diagnosis
        generalized_unifrac_fig_location_dir = os.path.join(
            beta_fig_location_dir, "generalized_unifrac"
        )

        # Path to NMDS plot file comparing locations per diagnosis
        nmds_plot_file_path = os.path.join(
            generalized_unifrac_fig_location_dir,
            f"generalized_unifrac_nmds_plot_comparing_locations_{diagnosis}.pdf",
        )

        if os.path.exists(nmds_plot_file_path):
            print(
                f"Generalized UniFrac NMDS plot already exists at {nmds_plot_file_path}. Skipping Generalized UniFrac analysis for {diagnosis}."
            )
        else:
            print(f"Calculating Generalized UniFrac distances for {diagnosis}...")

            # Path to Generalized UniFrac distance matrix comparing locations per diagnosis
            generalized_unifrac_distance_matrix_file = os.path.join(
                generalized_unifrac_location_dir,
                f"generalized_unifrac_distance_matrix_comparing_locations_{diagnosis}.csv",
            )

            # Calculate Generalized UniFrac distance
            generalized_unifrac_df, generalized_unifrac_dm = (
                calculate_generalized_unifrac_distance(
                    biom_table_path=biom_file_path,
                    rooted_tree_newick_path=rooted_tree_newick_path,
                    output_file=generalized_unifrac_distance_matrix_file,
                )
            )

            # Perform NMDS
            coordinates_df, normalized_stress = perform_nmds(
                distance_matrix=generalized_unifrac_dm
            )

            # Perform pairwise comparisons with vegan on the Generalized UniFrac distance matrix
            perform_pairwise_comparisons_permanova(
                metadata=diagnosis_df,
                distance_df=generalized_unifrac_df,
                primary_variable="sample_location",
                output_file=os.path.join(
                    generalized_unifrac_location_dir,
                    f"generalized_unifrac_pairwise_comparisons_between_locations_{diagnosis}.csv",
                ),
            )

            # Perform permanova  with vegan on the Generalized UniFrac distance matrix
            results_dict = perform_permanova_with_vegan(
                metadata=diagnosis_df,
                distance_df=generalized_unifrac_df,
                primary_variable="sample_location",
            )

            # Create NMDS plot
            title_dict = {
                "primary_variable": "sample_location",
                "title": f"NMDS Plot of Generalized UniFrac distances comparing locations in {diagnosis} Patients",
                "stress": normalized_stress,
            }

            # Plot NMDS
            plot_nmds_with_ggplot2(
                nmds_coordinates=coordinates_df,
                metadata=diagnosis_df,
                permanova_results=results_dict,
                title_dict=title_dict,
                output_file=nmds_plot_file_path,
            )

        ##############################################################
        # 3.2.4 Differential abundance between locations per diagnosis
        ##############################################################

        # Define directory for differential abundance results
        daa_dir = os.path.join(wor_dir, "differential_abundance_analysis")

        # Path to DAA results directory comparing locations per diagnosis
        location_daa_dir = os.path.join(daa_dir, "location")

        # Path to DAA figures directory
        daa_fig_dir = os.path.join(fig_dir, "differential_abundance_analysis")

        # Path to DAA figures directory for comparing locations per diagnosis
        location_daa_fig_dir = os.path.join(daa_fig_dir, "location")

        ####################
        # ASV level analysis
        ####################

        # Path to DAA figures directory for ASV level analysis comparing locations per diagnosis
        diagnosis_daa_asv_fig_dir = os.path.join(location_daa_fig_dir, "asv_level")

        # Figure path to save the DAA plot
        daa_asv_plot_file = os.path.join(
            diagnosis_daa_asv_fig_dir, f"daa_asv_heatmap_diagnosis_{diagnosis}.pdf"
        )

        if os.path.exists(daa_asv_plot_file):
            print(
                f"DAA at the ASV level for comparing locations for {diagnosis} already done. Skipping ASV level analysis."
            )
        else:
            print(f"Performing DAA at the ASV level for {diagnosis}...")

            # Path to DAA results directory for ASV level analysis comparing locations per diagnosis
            location_daa_asv_level = os.path.join(location_daa_dir, "asv_level")

            daa_file = os.path.join(
                location_daa_asv_level,
                f"differential_abundance_results_location_{diagnosis}.csv",
            )

            # Run LINDA analysis
            daa_results_taxa_fil_df = linda_daa_asv(
                asv_table=asv_table_diagnosis,
                metadata=diagnosis_df,
                primary_variable="sample_location",
                taxonomy_table=taxonomy_df,
                output_file=daa_file,
                threshold=5,
                pairwise=True,
            )

            # Plot the DAA results
            print("Plotting DAA results for ASVs...")

            # constitute title dictionary for the plot
            title_dict = {
                "primary_variable": "sample_location",
                "title": f"DAA Results for ASVs in {diagnosis} Patients",
            }

            plot_heatmap_of_daa(
                daa_df=daa_results_taxa_fil_df,
                title_dict=title_dict,
                output_file=daa_asv_plot_file,
                asv=True,  # Set asv to True for ASV level analysis
            )

        ########################
        # Taxa level analysis
        ########################

        # Taxonomic ranks to process
        ranks = ["species", "genus", "family", "order", "class", "phylum"]

        for rank in ranks:
            print(f"Processing taxonomic rank: {rank}")

            rank_daa_plot_file = os.path.join(
                location_daa_fig_dir,
                f"{rank}_level",
                f"daa_{rank}_heatmap_location_{diagnosis}.pdf",
            )

            if os.path.exists(rank_daa_plot_file):
                print(
                    f"DAA at the {rank} level for comparing locations for {diagnosis} already done. Skipping {rank} level analysis."
                )
            else:
                location_daa_rank_level = os.path.join(
                    location_daa_dir, f"{rank}_level"
                )

                # Define output file for DAA results
                daa_rank_file = os.path.join(
                    location_daa_rank_level,
                    f"daa_{rank}_results_location_{diagnosis}.csv",
                )

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

                taxon_daa_results_fil_df = linda_daa_taxon(
                    taxon_table_df=taxa_table_df,
                    metadata_df=diagnosis_df,
                    primary_variable="sample_location",
                    pairwise=True,
                    threshold=5,
                    output_file=daa_rank_file,
                )

                # Plot the DAA results
                print(f"Plotting DAA results for {rank} level...")

                # constitute title dictionary for the plot
                title_dict = {
                    "primary_variable": "sample_location",
                    "title": f"DAA between Locations at {rank} level for {diagnosis} Patients",
                }

                plot_heatmap_of_daa(
                    daa_df=taxon_daa_results_fil_df,
                    title_dict=title_dict,
                    output_file=rank_daa_plot_file,
                    asv=False,  # Set asv to False for taxonomic level analysis
                )

        print(f"Finished all analyses for diagnosis: {diagnosis}")

    print("Finished all diagnoses analyses.")

    ########################################################
    ########################################################
    # 3.3 Comparative analysis between diagnoses per location
    ########################################################
    ########################################################

    locations = merged_df["sample_location"].unique()

    for location in locations:
        print(f"Processing location: {location}")

        # Filter the merged DataFrame for the current location
        location_df = merged_df[merged_df["sample_location"] == location]

        # Filter the ASV table for the current location
        asv_table_location = asv_table_dep_fil_df[location_df.index]

        # Rarefy the ASV table
        rarefied_df = rarefy_df(asv_table_location)

        # Define output file for the ASV table in BIOM format
        biom_file_path = os.path.join(asv_table_dir, f"asv_table_{location}.biom")

        # Convert to BIOM format
        convert_to_biom_format(
            asv_table_df=asv_table_location,
            output_file=biom_file_path,
        )
        print(f"BIOM file saved to {biom_file_path}")

        ######################
        ######################
        # Ecological Diversity
        ######################
        ######################

        ######################
        # Alpha Diversity
        ######################

        # Path to results for figures directory comparing diagnoses per location
        diagnosis_alpha_fig_dir = os.path.join(alpha_fig_dir, "diagnosis")

        # Path to file to save alpha diversity figures
        alpha_diversity_figure_file = os.path.join(
            diagnosis_alpha_fig_dir,
            f"alpha_diversity_metrics_diagnosis_comparison_{location}.pdf",
        )

        if os.path.exists(alpha_diversity_figure_file):
            print(
                f"Alpha diversity figure already exists at {alpha_diversity_figure_file}. Skipping alpha diversity analysis for {location}."
            )
        else:
            print(f"Calculating alpha diversity for {location}...")

            # Subset the alpha diversity DataFrame for the current location
            alpha_diversity_location_df = alpha_diversity_df[
                alpha_diversity_df["sample_id"].isin(location_df.index)
            ]

            # Path to results of statiscal comparisons of alpha diversity between diagnoses per location
            alpha_diversity_stats_results = os.path.join(
                alpha_div_dir,
                "results",
                "diagnosis",
                f"alpha_diversity_diagnosis_comparison_{location}.csv",
            )

            # Calculate differences in alpha diversity between diagnoses
            alpha_diversity_results_df = unpaired_alpha_diversity_calculations(
                alpha_diversity_df=alpha_diversity_location_df,
                metadata_df=location_df,
                primary_variable="Diagnosis",
                output_file=alpha_diversity_stats_results,
            )

            plot_heatmap_of_alpha_diversity(
                data_df=alpha_diversity_results_df,
                title_dict={
                    "primary_variable": "Diagnosis",
                    "title": f"Alpha Diversity Comparison between Diagnoses for {location}",
                },
                output_file=alpha_diversity_figure_file,
            )

        #####################
        # Beta Diversity
        #####################

        # Path to beta diversity calculation results comparing diagnoses per location
        beta_div_diagnosis_dir = os.path.join(beta_div_dir, "diagnosis")

        # Path to beta diversity figures comparing diagnoses per location
        beta_fig_diagnosis_dir = os.path.join(beta_fig_dir, "diagnosis")

        #############################
        # 3.3.3.1 Jaccard
        #############################

        # Path to Jaccard figures comparing diagnoses per location
        jaccard_diagnosis_fig_dir = os.path.join(beta_fig_diagnosis_dir, "jaccard")

        # Path to NMDS plot file comparing locations per diagnosis
        nmds_plot_file_path = os.path.join(
            jaccard_diagnosis_fig_dir,
            f"jaccard_nmds_plot_comparing_diagnoses_{location}.pdf",
        )

        if os.path.exists(nmds_plot_file_path):
            print(
                f"Jaccard NMDS plot already exists at {nmds_plot_file_path}. Skipping Jaccard analysis for {location}."
            )
        else:
            print(f"Calculating Jaccard distances for {location}...")

            # Path to Jaccard directory
            jaccard_diagnosis_dir = os.path.join(beta_div_diagnosis_dir, "jaccard")

            # Path to Jaccard distance matrix comparing locations per diagnosis
            jaccard_distance_matrix_file = os.path.join(
                jaccard_diagnosis_dir,
                f"jaccard_distance_matrix_comparing_diagnoses_{location}.csv",
            )

            # Calculate Jaccard distance
            jaccard_df, jaccard_dm = calculate_jaccard_distance(
                rarefied_table_df=rarefied_df,
                output_file=jaccard_distance_matrix_file,
            )

            # Perform pairwise comparisons with vegan on the Jaccard distance matrix
            perform_pairwise_comparisons_permanova(
                metadata=location_df,
                distance_df=jaccard_df,
                primary_variable="Diagnosis",
                output_file=os.path.join(
                    jaccard_diagnosis_dir,
                    f"jaccard_pairwise_comparisons_between_diagnoses_{location}.csv",
                ),
            )

            # Perform permanova  with vegan on the Jaccard distance matrix comparing diagnoses per location
            permanova_results_dict = perform_permanova_with_vegan(
                metadata=location_df,
                distance_df=jaccard_df,
                primary_variable="Diagnosis",
            )

            # Perform NMDS on the Jaccard distance matrix
            coordinates_df, normalized_stress = perform_nmds(distance_matrix=jaccard_dm)

            # Define a title dictionary for the NMDS plot
            title_dict = {
                "primary_variable": "Diagnosis",
                "title": f"NMDS Plot of Jaccard distances comparing Diagnoses for {location}",
                "stress": normalized_stress,
            }

            plot_nmds_with_ggplot2(
                nmds_coordinates=coordinates_df,
                metadata=diagnosis_df,
                permanova_results=permanova_results_dict,
                title_dict=title_dict,
                output_file=nmds_plot_file_path,
            )

        #####################
        # 3.3.3.2 Bray-Curtis
        #####################

        # Path to Bray-Curtis figures comparing diagnoses per location
        beta_fig_diagnosis_bray_curtis_dir = os.path.join(
            beta_fig_diagnosis_dir, "bray_curtis"
        )

        # Define output file for NMDS plot
        nmds_plot_file_path = os.path.join(
            beta_fig_diagnosis_bray_curtis_dir,
            f"bray_curtis_nmds_plot_comparing_diagnoses_{location}.pdf",
        )

        if os.path.exists(nmds_plot_file_path):
            print(
                f"Bray-Curtis NMDS plot already exists at {nmds_plot_file_path}. Skipping Bray-Curtis analysis for {location}."
            )
        else:
            # Path to directory to results of Bray-Curtis comparing diagnoses per location
            bray_curtis_diagnosis_dir = os.path.join(
                beta_div_diagnosis_dir, "bray_curtis"
            )

            # Path to output file for Bray-Curtis distance matrix
            bray_curtis_dist_file = os.path.join(
                bray_curtis_diagnosis_dir,
                f"bray_curtis_distance_matrix_comparing_diagnoses_{location}.csv",
            )

            # Perform Bray-Curtis distance calculation
            bray_curtis_df, bray_curtis_dm = calculate_bray_curtis_distance(
                rarefied_table_df=rarefied_df, output_file=bray_curtis_dist_file
            )

            results_dict = perform_permanova_with_vegan(
                metadata=location_df,
                distance_df=bray_curtis_df,
                primary_variable="Diagnosis",
            )

            # Perform pairwise comparisons with vegan on the Bray-Curtis distance matrix
            perform_pairwise_comparisons_permanova(
                metadata=location_df,
                distance_df=bray_curtis_df,
                primary_variable="Diagnosis",
                output_file=bray_curtis_dist_file,
            )

            # Perform NMDS on the Bray-Curtis distance matrix
            coordinates_df, normalized_stress = perform_nmds(
                distance_matrix=bray_curtis_dm
            )

            # Define a title dictionary for the NMDS plot
            title_dict = {
                "primary_variable": "Diagnosis",
                "title": f"NMDS Plot of Bray-Curtis distances comparing Diagnoses for {location}",
                "stress": normalized_stress,
            }

            plot_nmds_with_ggplot2(
                nmds_coordinates=coordinates_df,
                metadata=location_df,
                title_dict=title_dict,
                permanova_results=results_dict,
                output_file=nmds_plot_file_path,
            )
            print(f"NMDS plot saved to {nmds_plot_file_path}")

        ##############################
        # 3.3.3.3 Robust Aitchison PCA
        ##############################

        # Path to figures directory for RPCA results comparing diagnoses per location
        rpca_fig_diagnosis_dir = os.path.join(beta_fig_diagnosis_dir, "rpca")

        rpca_biplot_file = os.path.join(
            rpca_fig_diagnosis_dir, f"rpca_biplot_comparing_diagnoses_{location}.pdf"
        )

        if os.path.exists(rpca_biplot_file):
            print(
                f"RPCA biplot already exists at {rpca_biplot_file}. Skipping RPCA analysis for {location}."
            )
        else:
            # Directory to results of Robust Aitchison PCA comparing diagnoses per location
            rpca_dir = os.path.join(beta_div_diagnosis_dir, "rpca")

            # define output file for RPCA distance matrix
            rpca_distance_matrix_file = os.path.join(
                rpca_dir, f"rpca_distance_matrix_comparing_diagnoses_{location}.csv"
            )

            # Perform RPCA
            rpca_ordination_results, robust_aitchison_distance_df = perform_rpca(
                asv_table_df=asv_table_location, output_file=rpca_distance_matrix_file
            )

            # Perform PERMANOVA on the RPCA distance matrix
            rpca_results_dict = perform_permanova_with_vegan(
                metadata=location_df,
                distance_df=robust_aitchison_distance_df,
                primary_variable="Diagnosis",
            )

            # Perform pairwise comparisons with vegan on the RPCA distance matrix
            perform_pairwise_comparisons_permanova(
                metadata=location_df,
                distance_df=robust_aitchison_distance_df,
                primary_variable="Diagnosis",
                output_file=os.path.join(
                    rpca_dir,
                    f"rpca_pairwise_comparisons_between_diagnoses_{location}.csv",
                ),
            )

            # Create a title dictionary for the plot
            title_dict = {
                "primary_variable": "Diagnosis",
                "title": f"RPCA Biplot for location {location}",
            }

            plot_biplot_with_ggplot2(
                ordination_results=rpca_ordination_results,
                metadata=location_df,
                title_dict=title_dict,
                permanova_results=rpca_results_dict,
                taxonomy_table=taxonomy_df,
                output_file=rpca_biplot_file,
            )

        ############
        # Phylo RPCA
        ############

        # Path to figures directory for Phylo RPCA results comparing diagnoses per location
        phylo_rpca_fig_diagnosis_dir = os.path.join(
            beta_fig_diagnosis_dir, "phylo_rpca"
        )

        biplot_file = os.path.join(
            phylo_rpca_fig_diagnosis_dir,
            f"phylo_rpca_biplot_comparing_diagnoses_{location}.pdf",
        )

        if os.path.exists(biplot_file):
            print(
                f"Phylogenetic RPCA biplot already exists at {biplot_file}. Skipping Phylogenetic RPCA analysis for {location}."
            )
        else:
            print(f"Calculating Phylogenetic RPCA for {location}...")

            # Path to directory to results of Phylogenetic Robust Aitchison PCA comparing diagnoses per location
            phylo_rpca_dir = os.path.join(
                beta_fig_diagnosis_dir, "phylogenetic_robust_aitchison_pca"
            )

            # Perform Phylogenetic RPCA
            (
                phylo_ordination,
                phylo_rpca_distance_df,
                pruned_tree,
                filtered_table,
                filtered_taxonomy,
            ) = perform_phylogenetic_rpca(
                asv_table_df=asv_table_location,
                rooted_tree_newick_path=rooted_tree_newick_path,
                taxonomy_df=taxonomy_df,
                output_file=os.path.join(
                    phylo_rpca_dir,
                    f"phylo_rpca_location_{location}_distance_matrix.csv",
                ),
            )

            # Perform PERMANOVA on the phylo RPCA distance matrix
            phylo_rpca_results_dict = perform_permanova_with_vegan(
                metadata=location_df,
                distance_df=phylo_rpca_distance_df,
                primary_variable="Diagnosis",
            )

            # Perform pairwise comparisons with vegan on the phylo RPCA distance matrix
            perform_pairwise_comparisons_permanova(
                metadata=location_df,
                distance_df=phylo_rpca_distance_df,
                primary_variable="Diagnosis",
                output_file=os.path.join(
                    phylo_rpca_dir,
                    f"phylo_rpca_pairwise_comparisons_between_diagnoses_{location}.csv",
                ),
            )

            # Create a title dictionary for the phylo RPCA plot
            phylo_title_dict = {
                "primary_variable": "Diagnosis",
                "title": f"Phylo RPCA Biplot comparing Diagnoses for location {location}",
            }

            plot_biplot_with_ggplot2(
                ordination_results=phylo_ordination,
                metadata=location_df,
                title_dict=phylo_title_dict,
                permanova_results=phylo_rpca_results_dict,
                taxonomy_table=filtered_taxonomy,
                output_file=biplot_file,
            )

        ##############################
        # 3.3.3.5  Unweighted UniFrac
        ##############################

        # Path to figures directory for Unweighted UniFrac results comparing locations per diagnosis
        unweighted_unifrac_fig_diagnosis_dir = os.path.join(
            beta_fig_diagnosis_dir, "unweighted_unifrac"
        )

        # Path to NMDS plot file comparing locations per diagnosis
        nmds_plot_file_path = os.path.join(
            unweighted_unifrac_fig_diagnosis_dir,
            f"unweighted_unifrac_nmds_plot_comparing_diagnoses_{location}.pdf",
        )

        if os.path.exists(nmds_plot_file_path):
            print(
                f"Unweighted UniFrac NMDS plot already exists at {nmds_plot_file_path}. Skipping Unweighted UniFrac analysis for {location}."
            )
        else:
            print(f"Calculating Unweighted UniFrac distances for {location}...")

            # Path to Unweighted UniFrac directory
            unweighted_unifrac_diagnosis_dir = os.path.join(
                beta_div_diagnosis_dir, "unweighted_unifrac"
            )

            # Path to Unweighted UniFrac distance matrix comparing locations per diagnosis
            unweighted_unifrac_distance_matrix_file = os.path.join(
                unweighted_unifrac_diagnosis_dir,
                f"unweighted_unifrac_distance_matrix_comparing_locations_{diagnosis}.csv",
            )

            # Calculate Unweighted UniFrac distance
            unweighted_unifrac_df, unweighted_unifrac_dm = (
                calculate_unweighted_unifrac_distance(
                    biom_table_path=biom_file_path,
                    rooted_tree_newick_path=rooted_tree_newick_path,
                    output_file=unweighted_unifrac_distance_matrix_file,
                )
            )

            # Perform NMDS
            coordinates_df, normalized_stress = perform_nmds(
                distance_matrix=unweighted_unifrac_dm
            )

            # Perform pairwise comparisons with vegan on the Unweighted UniFrac distance matrix
            perform_pairwise_comparisons_permanova(
                metadata=location_df,
                distance_df=unweighted_unifrac_df,
                primary_variable="Diagnosis",
                output_file=os.path.join(
                    unweighted_unifrac_diagnosis_dir,
                    f"unweighted_unifrac_pairwise_comparisons_between_diagnoses_{location}.csv",
                ),
            )

            # Perform permanova  with vegan on the Unweighted UniFrac distance matrix
            results_dict = perform_permanova_with_vegan(
                metadata=location_df,
                distance_df=unweighted_unifrac_df,
                primary_variable="Diagnosis",
            )

            # Create NMDS plot
            title_dict = {
                "primary_variable": "Diagnosis",
                "title": f"NMDS Plot of Unweighted UniFrac distances comparing diagnoses for {location}",
                "stress": normalized_stress,
            }

            # Plot NMDS
            plot_nmds_with_ggplot2(
                nmds_coordinates=coordinates_df,
                metadata=location_df,
                permanova_results=results_dict,
                title_dict=title_dict,
                output_file=nmds_plot_file_path,
            )

        ###############################
        # 3.3.3.5  Generalized UniFrac
        ##############################

        # Path to figures directory for Generalized UniFrac results comparing locations per diagnosis
        generalized_unifrac_fig_diagnosis_dir = os.path.join(
            beta_fig_diagnosis_dir, "generalized_unifrac"
        )

        # Path to NMDS plot file comparing locations per diagnosis
        nmds_plot_file_path = os.path.join(
            generalized_unifrac_fig_diagnosis_dir,
            f"generalized_unifrac_nmds_plot_comparing_diagnoses_{location}.pdf",
        )

        if os.path.exists(nmds_plot_file_path):
            print(
                f"Generalized UniFrac NMDS plot already exists at {nmds_plot_file_path}. Skipping Generalized UniFrac analysis for {location}."
            )
        else:
            print(f"Calculating Generalized UniFrac distances for {location}...")

            # Path to Generalized UniFrac directory
            generalized_unifrac_diagnosis_dir = os.path.join(
                beta_div_diagnosis_dir, "generalized_unifrac"
            )

            # Path to Generalized UniFrac distance matrix comparing locations per diagnosis
            generalized_unifrac_distance_matrix_file = os.path.join(
                generalized_unifrac_diagnosis_dir,
                f"generalized_unifrac_distance_matrix_comparing_diagnoses_{location}.csv",
            )

            # Calculate Generalized UniFrac distance
            generalized_unifrac_df, generalized_unifrac_dm = (
                calculate_generalized_unifrac_distance(
                    biom_table_path=biom_file_path,
                    rooted_tree_newick_path=rooted_tree_newick_path,
                    output_file=generalized_unifrac_distance_matrix_file,
                )
            )

            # Perform NMDS
            coordinates_df, normalized_stress = perform_nmds(
                distance_matrix=generalized_unifrac_dm
            )

            pairwise_output_file = os.path.join(
                generalized_unifrac_diagnosis_dir,
                f"generalized_unifrac_pairwise_comparisons_between_diagnoses_{location}.csv",
            )

            # Perform pairwise comparisons with vegan on the Generalized UniFrac distance matrix
            perform_pairwise_comparisons_permanova(
                metadata=location_df,
                distance_df=generalized_unifrac_df,
                primary_variable="Diagnosis",
                output_file=pairwise_output_file,
            )

            # Perform permanova  with vegan on the Generalized UniFrac distance matrix
            results_dict = perform_permanova_with_vegan(
                metadata=location_df,
                distance_df=generalized_unifrac_df,
                primary_variable="Diagnosis",
            )

            # Create NMDS plot
            title_dict = {
                "primary_variable": "Diagnosis",
                "title": f"NMDS Plot of Generalized UniFrac distances comparing locations in {diagnosis} Patients",
                "stress": normalized_stress,
            }

            # Plot NMDS
            plot_nmds_with_ggplot2(
                nmds_coordinates=coordinates_df,
                metadata=diagnosis_df,
                permanova_results=results_dict,
                title_dict=title_dict,
                output_file=nmds_plot_file_path,
            )

        ##############################
        # Machine Learning Analysis
        ###############################

        # Make new metadata column by grouping Diagnosis variables
        location_df["Diagnosis_Group"] = location_df["Diagnosis"].replace(
            {
                "Healthy": "pre-transformation",
                "GORD": "pre-transformation",
                "BO": "pre-transformation",
                "Dysplasia": "post-transformation",
                "OAC": "post-transformation",
                "Metastatic": "post-transformation",
            }
        )

        # make a 0-1 column for Diagnosis_Group for machine learning
        location_df["Diagnosis_Group_0_1"] = location_df["Diagnosis_Group"].replace(
            {
                "pre-transformation": 0,
                "post-transformation": 1,
            }
        )

        # model, model_results_dict = train_xgboost_model(
        #     count_table_df=asv_table_location,
        #     metadata_df=location_df,
        #     primary_label="Diagnosis_Group_0_1",
        #     normalization_method="presence_absence",
        # )

        print("Model training completed.")

        #########################################
        #########################################
        # Differential Abundance Analysis (LINDA)
        #########################################
        #########################################

        # Path to DAA results directory comparing locations per diagnosis
        daa_diagnosis_dir = os.path.join(daa_dir, "diagnosis")

        # Path to DAA results directory for ASV level analysis comparing locations per diagnosis

        daa_diagnosis_asv_level = os.path.join(daa_diagnosis_dir, "asv_level")

        daa_file = os.path.join(
            daa_diagnosis_asv_level,
            f"differential_abundance_results_comparing_diagnoses_{location}.csv",
        )

        # Path to DAA figures directory
        daa_fig_dir = os.path.join(fig_dir, "differential_abundance_analysis")

        # Path to DAA figures directory for comparing locations per diagnosis
        diagnoses_daa_fig_dir = os.path.join(daa_fig_dir, "diagnosis")

        # Path to DAA figures directory for ASV level analysis comparing locations per diagnosis
        diagnosis_daa_asv_fig_dir = os.path.join(diagnoses_daa_fig_dir, "asv_level")

        # Figure path to save the DAA plot
        daa_asv_plot_file = os.path.join(
            diagnosis_daa_asv_fig_dir,
            f"daa_asv_heatmap_comparing_diagnoses_{location}.pdf",
        )

        # Run LINDA analysis
        daa_results_taxa_fil_df = linda_daa_asv(
            asv_table=asv_table_location,
            metadata=location_df,
            primary_variable="Diagnosis",
            taxonomy_table=taxonomy_df,
            output_file=daa_file,
            threshold=5,
            pairwise=False,
        )

        # Plot the DAA results
        print("Plotting DAA results for ASVs...")

        # constitute title dictionary for the plot
        title_dict = {
            "primary_variable": "Diagnosis",
            "title": f"Heatmap of DAA Results for ASVs comparing Diagnoses for location {location}",
        }

        plot_heatmap_of_daa(
            daa_df=daa_results_taxa_fil_df,
            title_dict=title_dict,
            output_file=daa_asv_plot_file,
            asv=True,
        )

        # Taxonomic ranks to process
        ranks = ["species", "genus", "family", "order", "class", "phylum"]

        for rank in ranks:
            print(f"Processing taxonomic rank: {rank}")

            rank_daa_plot_file = os.path.join(
                wor_dir,
                "figures",
                "differential_abundance_analysis",
                "diagnosis",
                f"{rank}_level",
                f"daa_{rank}_heatmap_comparing_diagnoses_{location}.pdf",
            )

            if os.path.exists(rank_daa_plot_file):
                print(
                    f"DAA {rank} plot already exists at for location {location}. Skipping DAA analysis for {rank} level for location {location}."
                )
            else:
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
                    "differential_abundance_analysis",
                    "diagnosis",
                    f"{rank}_level",
                    f"daa_{rank}_results_comparing_diagnoses_{location}.csv",
                )

                taxon_daa_results_fil_df = linda_daa_taxon(
                    taxa_table_df,
                    location_df,
                    primary_variable="Diagnosis",
                    pairwise=False,
                    threshold=5,
                    output_file=output_file,
                )

                # Plot the DAA results
                print(f"Plotting DAA results for {rank} level...")

                # constitute title dictionary for the plot
                title_dict = {
                    "primary_variable": "Diagnosis",
                    "title": f"DAA Results for {rank} level at location {location}",
                }

                plot_heatmap_of_daa(
                    daa_df=taxon_daa_results_fil_df,
                    title_dict=title_dict,
                    output_file=rank_daa_plot_file,
                    asv=False,  # Set asv to False for taxonomic level analysis
                )

                print(f"Finished processing {rank} level for location {location}.")

        ##########################################################################
        ##########################################################################
        # Differential Abundance Analysis (LINDA) on predicted functional features
        ##########################################################################
        ##########################################################################

        # Path to PICRUSt2 analysis directory
        picrust2_dir = os.path.join(wor_dir, "picrust2_differential_abundance_analysis")

        # Path to PICRUSt2 analysis results comparing diagnoses per location
        picrust2_diagnosis_dir = os.path.join(picrust2_dir, "diagnosis")

        # Path to figures for Picrust2 analysis
        picrust2_fig_dir = os.path.join(
            fig_dir, "picrust2_differential_abundance_analysis"
        )

        # Path to figures for Picrust2 analysis comparing diagnoses per location
        picrust2_diagnosis_fig_dir = os.path.join(picrust2_fig_dir, "diagnosis")

        print("Performing DAA on pathways...")

        # Define all the pathway types and their corresponding variables
        pathway_info = [
            {
                "name": "KEGG",
                "table_df": kegg_pathways_df,
                "subfolder": "kegg_level",
                "file_prefix": "kegg",
            },
            {
                "name": "EC",
                "table_df": enzyme_commission_df,
                "subfolder": "ec_level",
                "file_prefix": "ec",
            },
            {
                "name": "MetaCyc",
                "table_df": metabolic_pathways_df,
                "subfolder": "metacyc_level",
                "file_prefix": "metacyc",
            },
        ]

        for info in pathway_info:
            print(f"Performing DAA on {info['name']} pathways...")

            table_location = info["table_df"][location_df.index]

            # Define output file for DAA results
            daa_file = os.path.join(
                wor_dir,
                picrust2_diagnosis_dir,
                info["subfolder"],
                f"daa_{info['file_prefix']}_results_comparing_diagnoses_{location}.csv",
            )

            # Run LINDA analysis
            daa_results_fil_df = linda_daa_picrust2(
                function_table_df=table_location,
                metadata_df=location_df,
                primary_variable="Diagnosis",
                threshold=20,
                pairwise=False,
                output_file=daa_file,
            )

            # Plot the DAA results
            print(f"Plotting DAA results for {info['name']} pathways...")
            plot_file_path = os.path.join(
                wor_dir,
                picrust2_diagnosis_fig_dir,
                info["subfolder"],
                f"daa_{info['file_prefix']}_heatmap_comparing_diagnoses_{location}.pdf",
            )

            # constitute title dictionary for the plot
            title_dict = {
                "primary_variable": "Diagnosis",
                "title": f"DAA Results for {info['name']} functions comparing Diagnoses for {location}",
            }

            plot_heatmap_of_daa(
                daa_df=daa_results_fil_df,
                title_dict=title_dict,
                output_file=plot_file_path,
                asv=False,  # Set asv to False for pathway level analysis
            )

        print(f"Finished all analyses for location: {location}")

    print("Finished all locations analyses.")

    ##################
    # Network analysis
    ##################

    # ranks = ["species", "genus", "family", "order", "class", "phylum"]

    # for rank in ranks:
    #     # Define output file path for taxa tables
    #     taxa_tables_file = os.path.join(
    #         wor_dir, "tables", "taxa_tables", f"{rank}_table.csv"
    #     )

    #     # Create taxa tables
    #     taxa_table_df = create_taxa_tables(
    #         asv_table_df,
    #         taxonomy_df,
    #         rank=rank,
    #         output_file=taxa_tables_file,
    #     )

    #     # Define title dictionary for co-occurrence analysis
    #     title_dict = {
    #         "primary_variable": "Diagnosis",
    #         "location": str(location),  # Convert to string here
    #         "rank": rank,
    #     }

    #     # Define output directory for co-occurrence analysis
    #     output_dir = os.path.join(
    #         wor_dir,
    #         "co_occurrence_analysis",
    #         f"{rank}_level",
    #         str(location),  # Convert to string here too
    #     )

    #     print("Performing co-occurrence analysis for rank:", rank)

    #     perform_co_occurrence_network(
    #         taxa_table_df=taxa_table_df,
    #         title_dict=title_dict,
    #         output_dir=output_dir,
    # )

    # print(f"Finished processing location: {location}")

    # End of location loop

    # Copy the results to Dropbox
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
