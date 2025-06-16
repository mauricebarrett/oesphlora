import os
import subprocess
import sys
from pathlib import Path

import bionumpy as bnp  # type: ignore
import matplotlib.pyplot as plt
import numpy as np


def count_nreads_fastq(fastq_file, threads, output_file):
    """
    Count the number of reads in a FASTQ file using fqkit size command.

    Parameters:
    fastq_file (str): Path to the FASTQ file.
    threads (int): Number of threads to use for counting.
    output_file (str): Path to save the read count result.
    """

    # Skip if output file already exists
    if os.path.exists(output_file):
        print(f"Output file {output_file} already exists. Skipping read count.")
        return

    # Check if the file exists
    if not os.path.isfile(fastq_file):
        print(f"Error: The file {fastq_file} does not exist.")
        sys.exit(1)

    # Run seqkit stats command to count reads
    cmd = ["fqkit", "size", "--threads", str(threads), fastq_file]
    try:
        # Create the output directory if it doesn't exist
        Path(os.path.dirname(output_file)).mkdir(parents=True, exist_ok=True)
        # Run the command and capture the output
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        # Save the read count to the output file
        with open(output_file, "w") as f:
            f.write(result.stdout)
        print(f"Read count saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(cmd)}")
        print(f"Error message: {e.stderr}")
        sys.exit(1)


def gather_read_length_distribution(fastq_file, threads, output_file):
    """
    Glean read length distribution from a FASTQ file using fqkit length

    Parameters:
    fastq_file (str): Path to the FASTQ file.
    threads (int): Number of threads to use for processing.
    output_file (str): Path to save the read length distribution result.
    """

    # Skip if output file already exists
    if os.path.exists(output_file):
        print(f"Output file {output_file} already exists. Skipping read count.")
        return

    # Check if the file exists
    if not os.path.isfile(fastq_file):
        print(f"Error: The file {fastq_file} does not exist.")
        sys.exit(1)

    # Run fqkit length command to get read lengths

    cmd = ["fqkit", "length", "--threads", str(threads), fastq_file]
    try:
        # Create the output directory if it doesn't exist
        Path(os.path.dirname(output_file)).mkdir(parents=True, exist_ok=True)
        # Run the command and capture the output
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

        with open(output_file, "w") as f:
            f.write(result.stdout)
        print(f"Read length distribution saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(cmd)}")
        print(f"Error message: {e.stderr}")
        sys.exit(1)


def graph_quality_scores(fastq_file, threads, output_file):
    """
    Graph quality scores from a FASTQ file

    Parameters:
    fastq_file (str): Path to the FASTQ file.
    threads (int): Number of threads to use for processing.
    output_file (str): Path to save the quality score graph result.
    """

    # Check if the output file already exists
    if os.path.exists(output_file):
        print(f"Output file {output_file} already exists. Skipping graph generation.")
        return

    # Check if the file exists
    if not os.path.isfile(fastq_file):
        print(f"Error: The file {fastq_file} does not exist.")
        sys.exit(1)

    # Define subset of fastq file to read
    if isinstance(fastq_file, Path):
        fastq_file_str = str(fastq_file)
    else:
        fastq_file_str = fastq_file
    subset_fastq_file = fastq_file_str.replace(".fastq.gz", "_subset.fastq")

    sample_name = os.path.basename(fastq_file_str).split(".")[0]

    # Check if the output directory exists, if not create it
    output_dir = os.path.dirname(output_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Extract 2000 reads from the fastq file
    subprocess.run(
        f"fqkit topn -n 1000 --threads {threads} {fastq_file_str} > {subset_fastq_file}",
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

    # Save the figure
    plt.savefig(output_file, bbox_inches="tight")
    plt.close()


def fastq_qc(fastq_file, threads, out_dir):
    """
    Run all analysis steps on a FASTQ file.

    Parameters:
    fastq_file (str or Path): Path to the FASTQ file.
    threads (int): Number of threads.
    out_dir (str or Path): Directory for output files.
    """

    # Ensure out_dir is a Path object
    out_dir = Path(out_dir)

    # Get the base name of the FASTQ file
    base_name = os.path.basename(fastq_file).replace(".fastq.gz", "")

    # Define output files
    nreads_file = out_dir / f"{base_name}_readcount.txt"
    lengthdist_file = out_dir / f"{base_name}_lengthdist.txt"
    qualplot_file = out_dir / f"{base_name}_qualplot.png"

    # Call your analysis functions
    count_nreads_fastq(str(fastq_file), threads, str(nreads_file))
    gather_read_length_distribution(str(fastq_file), threads, str(lengthdist_file))
    graph_quality_scores(str(fastq_file), threads, str(qualplot_file))
