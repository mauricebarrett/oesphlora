if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install GenomicInfoDbData if not installed
if (!requireNamespace("GenomeInfoDbData", quietly = TRUE)) {
  BiocManager::install("GenomeInfoDbData")
}


a # Load required libraries
library(optparse)
library(dada2)
# Define command-line options
option_list <- list(
  make_option(c("-w", "--working_directory"),
    type = "character", default = NULL,
    help = "Working directory path", metavar = "character"
  )
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate input arguments
if (is.null(opt$working_dir)) {
  print_help(opt_parser)
  stop("--working_directory argument must be provided.",
    call. = FALSE
  )
}

# Assign variables
wor_dir <- opt$working_directory

# Main pipeline

# Define path to fastq files
fastq_files_path <- file.path(
  wor_dir, "fastq_files",
  "primers_removed"
)

# List forward and reverse FASTQ files
fwd_fqs <- sort(list.files(fastq_files_path,
  pattern = "_R1.fastq.gz",
  full.names = TRUE
))
rev_fqs <- sort(list.files(fastq_files_path,
  pattern = "_R2.fastq.gz",
  full.names = TRUE
))

# Extract sample names
sample_names <- gsub("_nop.*", "", basename(fwd_fqs))

# Filtered file paths
fwd_filt_path <- file.path(
  wor_dir, "fastq_files", "filtered",
  paste0(sample_names, "_F_filt.fastq.gz")
)
rev_filt_path <- file.path(
  wor_dir, "fastq_files", "filtered",
  paste0(sample_names, "_R_filt.fastq.gz")
)
names(fwd_filt_path) <- sample_names
names(rev_filt_path) <- sample_names

# Filter and trim reads
out <- filterAndTrim(fwd_fqs, fwd_filt_path, rev_fqs, rev_filt_path,
  truncLen = c(260, 230), maxN = 0, maxEE = c(1, 1), truncQ = 2,
  rm.phix = TRUE, compress = TRUE, multithread = FALSE
)

# Learn error rates
err_fwd <- learnErrors(fwd_filt_path, multithread = TRUE, verbose = TRUE)
err_rev <- learnErrors(rev_filt_path, multithread = TRUE, verbose = TRUE)

# Sample inference and merging
dada_fwd <- dada(fwd_filt_path, err = err_fwd, multithread = TRUE)
dada_rev <- dada(rev_filt_path, err = err_rev, multithread = TRUE)

mergers <- mergePairs(dada_fwd, fwd_filt_path,
  dada_rev, rev_filt_path,
  verbose = TRUE
)

# Create sequence table
seqtab <- makeSequenceTable(mergers)

# Remove sequences less than 350 and longer than 470
seqtab2 <- seqtab[, nchar(colnames(seqtab)) %in% 350:470]

# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab2,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

# Track reads through the pipeline
get_num <- function(x) sum(getUniques(x))
track <- cbind(
  out,
  sapply(dada_fwd, get_num),
  sapply(dada_rev, get_num),
  sapply(mergers, get_num),
  rowSums(seqtab),
  rowSums(seqtab2),
  rowSums(seqtab_nochim)
)
colnames(track) <- c(
  "input", "filtered", "denoisedF",
  "denoisedR", "merged", "seqtab", "seqtab2", "nonchim"
)

# Declare QC directory
qc_dir <- file.path(wor_dir, "qc", "read_count")

# Write out read count figures
write.table(track,
  file.path(qc_dir, "read_counts.csv"),
  sep = ",", quote = FALSE, row.names = TRUE
)

# Rename columns to sample names
rownames(track) <- sample_names
head(track, n = 20)

# Remove samples with read number below 5000
valid_samples <- names(which(rowSums(seqtab_nochim) > 4000))
seqtab_nochim_5000 <- as.matrix(seqtab_nochim[valid_samples, ])

# Write representative sequences to FASTA
sequences <- getSequences(seqtab_nochim_5000)

# Determine the fixed length based on the number of ASVs
num_asvs <- ncol(seqtab_nochim_5000)
fixed_length <- nchar(as.character(num_asvs))

# Generate column names with leading zeros
sequence_labels <- sprintf(
  paste0(">asv_%0", fixed_length, "d"),
  seq_len(num_asvs)
)

fasta_content <- paste(sequence_labels, sequences, sep = "\n")
output_dir <- file.path(wor_dir, "representative_sequences")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Create fasta file path
fasta_file <- file.path(
  output_dir,
  "asv_sequences.fasta"
)

writeLines(fasta_content, fasta_file)
cat("FASTA file written to:", fasta_file, "\n")

# Generate column names with leading zeros
colnames_labels <- sprintf(
  paste0("asv_%0", fixed_length, "d"),
  seq_len(num_asvs)
)

asv_table <- seqtab_nochim_5000
colnames(asv_table) <- colnames_labels

output_dir <- file.path(wor_dir, "tables", "asv_tables")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Write out ASV table
write.table(t(asv_table),
  file.path(output_dir, "asv_table.csv"),
  sep = ",", quote = FALSE, row.names = TRUE
)
