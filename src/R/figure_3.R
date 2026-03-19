library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)


# Set working directory (match analysis output location)
work_dir <- "/home/mossy/projects/phd/oesphlora"

# Paths to DAA results (comparing diagnoses per biopsy location)
daa_diagnosis_dir <- file.path(
  work_dir,
  "differential_abundance_analysis",
  "diagnosis"
)

# Biopsy site of interest for the main Figure 2
location <- "biopsy_location_3"

# ASV-level filtered DAA results for this site
asv_results_file <- file.path(
  daa_diagnosis_dir,
  "asv_level",
  paste0(
    "differential_abundance_results_comparing_diagnoses_",
    location,
    "_filtered.csv"
  )
)

if (!file.exists(asv_results_file)) {
  stop(
    paste(
      "ASV-level DAA results file not found for",
      location,
      "at:",
      asv_results_file
    )
  )
}

daa_df <- read.csv(asv_results_file, stringsAsFactors = FALSE)

# If taxonomy columns are present, create a more informative feature label
taxonomic_ranks <- c("phylum", "class", "order", "family", "genus", "species")

if (all(taxonomic_ranks %in% colnames(daa_df))) {
  lowest_classified <- apply(
    daa_df[, taxonomic_ranks],
    1,
    function(row) {
      val <- NA
      for (v in rev(row)) {
        if (!is.na(v) && v != "unclassified") {
          val <- v
          break
        }
      }
      if (is.na(val)) {
        "unclassified"
      } else {
        val
      }
    }
  )

  daa_df$feature_id <- paste0(daa_df$feature_id, "=", lowest_classified)
}

# Title information (matches usage in daa_heatmap.R); no plot title for paper figures
title_list <- list(
  "Diagnosis",
  ""
)

# Output path for main Figure 2
output_dir <- file.path(work_dir, "figures", "main_figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_file <- file.path(output_dir, "figure_3.pdf")

# Build heatmap directly (standalone, no external source files)
daa_df$comparison <- gsub("biopsy_location_", "", daa_df$comparison)
group_levels <- c("Healthy", "GORD", "BO", "Dysplasia", "OAC", "Metastatic")
all_comparisons <- combn(group_levels, 2, function(x) paste(x[1], "vs", x[2]))

all_combos <- expand.grid(
  feature_id = unique(daa_df$feature_id),
  comparison = all_comparisons,
  stringsAsFactors = FALSE
)

plot_full_df <- merge(
  daa_df,
  all_combos,
  by = c("comparison", "feature_id"),
  all = TRUE
)

plot_full_df$log2FoldChange[is.na(plot_full_df$log2FoldChange)] <- 0
plot_full_df$padj[is.na(plot_full_df$padj)] <- 1

plot_full_df <- plot_full_df %>%
  mutate(
    star = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  select(feature_id, comparison, log2FoldChange, star, padj)

data_matrix <- dcast(
  plot_full_df,
  feature_id ~ comparison,
  value.var = "log2FoldChange"
)
rownames(data_matrix) <- data_matrix$feature_id
data_matrix <- data_matrix[, -1]
data_matrix[is.na(data_matrix)] <- 0

if (nrow(data_matrix) > 1) {
  hc <- hclust(dist(data_matrix, method = "euclidean"), method = "ward.D2")
  ordered_features <- rownames(data_matrix)[hc$order]
  plot_full_df$feature_id <- factor(plot_full_df$feature_id, levels = ordered_features)
}

plot_full_df$comparison_standard <- factor(
  plot_full_df$comparison,
  levels = all_comparisons
)

figure_2 <- ggplot(plot_full_df, aes(x = comparison_standard, y = feature_id)) +
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
    size = 4,
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
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(colour = "black", size = 10, face = "bold"),
    legend.text = element_text(colour = "black", size = 10, face = "italic"),
    plot.caption = element_text(face = "bold")
  )

print(figure_2)

pdf(output_file, width = 14, height = 7.88, useDingbats = FALSE)
print(figure_2)
dev.off()

