#!/bin/bash

# Retrive qiiime2-amplicon-2024.10-py310-linux-conda.yml
wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml

# Convert yml to pixi toml
pixi init --import qiime2-amplicon-2024.10-py310-linux-conda.yml

# Install the environment
pixi install

# Add a package to the environment
pixi add --pypi "bionumpy==1.0.13"

pixi add cutadapt==5.0

pixi add chopper==0.9.2

pixi add seqkit==2.10.0

pixi add polars=1.26.0

pixi add "slack-sdk=3.35.0"

pixi add --pypi pyarrow

pixi add r-ggpubr

# Add a channel to the environment
pixi project channel add bioconda

pixi project channel add bioconda conda-forge

# Set the source and destination directories
SOURCE_DIR="/home/maurice/projects/phd/oesphlora/data/origin/oac_biopsies"
DEST_DIR="/home/maurice/projects/phd/oesphlora/fastq_files/demultiplexed"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Loop through run directories in sorted order
for run_dir in $(ls -d "$SOURCE_DIR"/run* 2>/dev/null | sort -V); do
    echo "Processing $run_dir"
    src_fastq_dir="$run_dir/fastq"
    if [ -d "$src_fastq_dir" ]; then
        cp "$src_fastq_dir"/*.fastq.gz "$DEST_DIR"/
    fi
done
