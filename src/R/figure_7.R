library(ggplot2)
library(dplyr)
library(reshape2)


# Match working directory to analysis output location
work_dir <- "/home/mossy/projects/phd/oesphlora"

# Diagnosis of interest for Figure 7 (Barrett's esophagus)
diagnosis <- "BO"

# Differential abundance ASV results comparing biopsy locations within diagnosis
daa_results_dir <- file.path(
  work_dir,
  "differential_abundance_analysis",
  "location",
  "asv_level"
)

filtered_results_file <- file.path(
  daa_results_dir,
  paste0("differential_abundance_results_location_", diagnosis, "_filtered.csv")
)

unfiltered_results_file <- file.path(
  daa_results_dir,
  paste0("differential_abundance_results_location_", diagnosis, ".csv")
)

if (file.exists(filtered_results_file)) {
  daa_results_file <- filtered_results_file
} else if (file.exists(unfiltered_results_file)) {
  daa_results_file <- unfiltered_results_file
} else {
  stop(
    paste(
      "Differential abundance ASV results file not found for diagnosis",
      diagnosis,
      "at:",
      filtered_results_file,
      "or",
      unfiltered_results_file
    )
  )
}

daa_df <- read.csv(daa_results_file, stringsAsFactors = FALSE)

if (nrow(daa_df) == 0) {
  stop("Differential abundance ASV results file is empty for diagnosis ", diagnosis, ".")
}

required_cols <- c("feature_id", "comparison", "log2FoldChange", "padj")
missing_cols <- setdiff(required_cols, colnames(daa_df))
if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required columns in differential abundance ASV results:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# If the script falls back to unfiltered results, retain significant ASVs only.
if (identical(daa_results_file, unfiltered_results_file)) {
  daa_df <- daa_df %>% dplyr::filter(!is.na(padj), padj < 0.05)
}

if (nrow(daa_df) == 0) {
  stop("No significant ASVs available to plot for diagnosis ", diagnosis, ".")
}

# Add lowest classified taxonomy label when columns are available.
taxonomy_ranks <- c("phylum", "class", "order", "family", "genus", "species")
if (all(taxonomy_ranks %in% colnames(daa_df))) {
  lowest_classified <- apply(
    daa_df[, taxonomy_ranks, drop = FALSE],
    1,
    function(row_values) {
      chosen <- NA_character_
      for (v in rev(row_values)) {
        if (!is.na(v) && nzchar(v) && v != "unclassified") {
          chosen <- as.character(v)
          break
        }
      }
      if (is.na(chosen)) {
        "unclassified"
      } else {
        chosen
      }
    }
  )
  daa_df$feature_label <- paste0(daa_df$feature_id, "=", lowest_classified)
} else {
  daa_df$feature_label <- daa_df$feature_id
}

# Create full pairwise-comparison grid so all cells are shown.
location_levels <- paste0("biopsy_location_", 1:5)
all_comparisons <- combn(location_levels, 2, function(x) paste(x[1], "vs", x[2]))

all_combos <- expand.grid(
  feature_label = unique(daa_df$feature_label),
  comparison = all_comparisons,
  stringsAsFactors = FALSE
)

plot_full_df <- merge(
  daa_df %>% dplyr::select(feature_label, comparison, log2FoldChange, padj),
  all_combos,
  by = c("comparison", "feature_label"),
  all = TRUE
)

plot_full_df$log2FoldChange[is.na(plot_full_df$log2FoldChange)] <- 0
plot_full_df$padj[is.na(plot_full_df$padj)] <- 1

plot_full_df <- plot_full_df %>%
  dplyr::mutate(
    star = dplyr::case_when(
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Drop biopsy_location_ prefix in axis labels
plot_full_df$comparison <- gsub("biopsy_location_", "", plot_full_df$comparison)
all_comparison_labels <- gsub("biopsy_location_", "", all_comparisons)

# Cluster ASVs by log2FC profile
data_matrix <- dcast(
  plot_full_df,
  feature_label ~ comparison,
  value.var = "log2FoldChange"
)
rownames(data_matrix) <- data_matrix$feature_label
data_matrix <- data_matrix[, -1, drop = FALSE]
data_matrix[is.na(data_matrix)] <- 0

if (nrow(data_matrix) > 1 && ncol(data_matrix) > 0) {
  hc <- hclust(dist(data_matrix, method = "euclidean"), method = "ward.D2")
  ordered_features <- rownames(data_matrix)[hc$order]
  plot_full_df$feature_label <- factor(plot_full_df$feature_label, levels = ordered_features)
}

plot_full_df$comparison_standard <- factor(
  plot_full_df$comparison,
  levels = all_comparison_labels
)

figure_7 <- ggplot(plot_full_df, aes(x = comparison_standard, y = feature_label)) +
  geom_tile(aes(fill = log2FoldChange), color = "black") +
  scale_fill_gradient2(
    low = "darkblue",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    space = "Lab",
    na.value = "grey90",
    name = "log2FC"
  ) +
  geom_text(
    aes(label = star),
    color = "black",
    size = 3,
    vjust = 0.7,
    na.rm = TRUE
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = NULL,
    caption = "* <= 0.05, ** <= 0.01, *** <= 0.001"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(colour = "black", size = 10, face = "bold"),
    legend.text = element_text(colour = "black", size = 10, face = "italic"),
    plot.caption = element_text(face = "bold"),
    plot.title = element_text(
      size = 14,
      hjust = 0.5,
      face = "bold",
      color = "black"
    )
  )

print(figure_7)

output_dir <- file.path(work_dir, "figures", "main_figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_file <- file.path(output_dir, "figure_7.pdf")

# Scale figure height with number of plotted ASVs for readability.
n_features <- length(unique(plot_full_df$feature_label))
figure_height <- max(8, min(18, 0.2 * n_features + 4))

pdf(output_file, width = 14, height = figure_height, useDingbats = FALSE)
print(figure_7)
dev.off()
