import argparse
import glob
import itertools
import os
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

import bionumpy as bnp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from biom.table import Table


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
def gather_qc_data(directory, stage):
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
    if fastq_files:
        print(f".fastq files found for {stage}, proceeding with processing.")
    else:
        print(f"No .fastq files found for {stage}. Exiting with error.")
        sys.exit(1)

    for filename in fastq_files:
        sample_name = os.path.splitext(filename)[0].split("_1")[0]
        print(f"Processing sample {sample_name}")

        # Counting number of reads in sample
        print(f"Counting number of reads in sample {sample_name}")
        fastq_file_path = os.path.join(fastq_dir, filename)
        count_cmd = ["seqkit", "stats", fastq_file_path]
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
            seqkit_cmd = ["seqkit", "fx2tab", "-nl", fastq_file_path]
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

    # Prepare a figure to save all the quality plots
    output_dir = os.path.join(directory, "qc", "q_score", stage)
    os.makedirs(output_dir, exist_ok=True)

    # Start generating quality plots
    for fastq_file in fastq_files:
        # Extract sample name
        sample_name = os.path.basename(fastq_file).split(".")[0]

        # Load FASTQ file with BioNumPy
        print(f"Graphing Phred Quality Score of {sample_name}, {stage} sample")
        f = bnp.open(fastq_file)
        data = f.read()  # This reads the sequences and quality scores

        # Extract the quality scores
        quality_scores = data.quality  # quality_scores is a NumPy array-like object

        # Sample a subset of reads if the dataset is too large
        num_reads_to_sample = 1000  # Adjust this value based on data
        quality_scores = (
            quality_scores[:num_reads_to_sample]
            if len(quality_scores) > num_reads_to_sample
            else quality_scores
        )

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
    directory, adapter_fwd="CCTACGGGNGGCWGCAG", adapter_rev="GACTACHVGGGTATCTAATCC"
):
    """
    Removes primers from paired-end FASTQ files using cutadapt.

    Parameters:
        directory (str): Root directory containing fastq_files/demultiplexed/
        adapter_fwd (str): Forward adapter sequence
        adapter_rev (str): Reverse adapter sequence
    """
    files_path = os.path.join(directory, "fastq_files", "demultiplexed")
    output_path = os.path.join(directory, "fastq_files", "trimmed")
    os.makedirs(output_path, exist_ok=True)

    r1_files = glob.glob(os.path.join(files_path, "*_R1_001.fastq.gz"))

    for r1 in r1_files:
        r2 = r1.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")

        if not os.path.exists(r2):
            print(f"⚠️  Missing R2 file for {os.path.basename(r1)}, skipping...")
            continue

        sample_name = os.path.basename(r1).replace("_R1_001.fastq.gz", "")
        out_r1 = os.path.join(output_path, f"{sample_name}_trimmed_R1.fastq.gz")
        out_r2 = os.path.join(output_path, f"{sample_name}_trimmed_R2.fastq.gz")

        cutadapt_cmd = [
            "cutadapt",
            "-a",
            adapter_fwd,
            "-A",
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


def main():
    # Parse the arguments using the parse_arguments function
    num_threads, directory = parse_arguments()

    # Keywords to filter
    keywords_to_filter = ["Chloroplast", "Mitochondrion", "Capsicum"]

    # Display the values for now (can be used later in the script)
    print(f"Number of threads: {num_threads}")
    print(f"Working directory: {directory}")

    # Gather QC data for each stage
    gather_qc_data(directory, "demultiplexed")

    # Graph quality scores of 16s FASTQ files
    graph_quality_scores(directory, "demultiplexed")


if __name__ == "__main__":
    main()
