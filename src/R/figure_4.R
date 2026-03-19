library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)


work_dir <- "/home/mossy/projects/phd/oesphlora"

picrust2_diagnosis_dir <- file.path(
  work_dir,
  "picrust2_differential_abundance_analysis",
  "diagnosis"
)

location <- "biopsy_location_3"

caption_text <- "* <= 0.05, ** <= 0.01, *** <= 0.001"
group_levels <- c("Healthy", "GORD", "BO", "Dysplasia", "OAC", "Metastatic")
all_comparisons <- combn(group_levels, 2, function(x) paste(x[1], "vs", x[2]))


prepare_daa_heatmap_df <- function(daa_df) {
  daa_df$comparison <- gsub("biopsy_location_", "", daa_df$comparison)

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

  if (nrow(data_matrix) > 1 && ncol(data_matrix) > 0) {
    hc <- hclust(dist(data_matrix, method = "euclidean"), method = "ward.D2")
    ordered_features <- rownames(data_matrix)[hc$order]
    plot_full_df$feature_id <- factor(
      plot_full_df$feature_id,
      levels = ordered_features
    )
  }

  plot_full_df$comparison <- factor(
    plot_full_df$comparison,
    levels = all_comparisons
  )

  plot_full_df
}


load_daa_data <- function(pathway_name, subfolder, file_prefix) {
  daa_file <- file.path(
    picrust2_diagnosis_dir,
    subfolder,
    paste0(
      "daa_",
      file_prefix,
      "_results_comparing_diagnoses_",
      location,
      "_filtered.csv"
    )
  )

  if (!file.exists(daa_file)) {
    stop(paste(
      "PICRUSt2 DAA results file not found for",
      pathway_name,
      "at:",
      daa_file
    ))
  }

  daa_df <- read.csv(daa_file, stringsAsFactors = FALSE)

  required_cols <- c("feature_id", "comparison", "log2FoldChange", "padj")
  if (nrow(daa_df) == 0 || !all(required_cols %in% colnames(daa_df))) {
    stop(paste("No plottable results for", pathway_name))
  }

  daa_df <- daa_df %>% filter(!is.na(feature_id), !is.na(comparison))
  if (nrow(daa_df) == 0 || length(unique(daa_df$feature_id)) < 2) {
    stop(paste("Not enough features to plot for", pathway_name))
  }

  prepare_daa_heatmap_df(daa_df)
}


ko_df <- load_daa_data(
  pathway_name = "KEGG KO",
  subfolder = "kegg_level",
  file_prefix = "kegg"
)

ec_df <- load_daa_data(
  pathway_name = "Enzyme Commission (EC)",
  subfolder = "ec_level",
  file_prefix = "ec"
)

facet_df <- rbind(
  cbind(ko_df, tag = "A"),
  cbind(ec_df, tag = "B")
)

panel_order <- c("A", "B")
facet_df$tag <- factor(facet_df$tag, levels = panel_order)
strip_labeller <- as_labeller(
  c(
    A = "KEGG Orthologues",
    B = "Enzyme Commission"
  )
)

tag_df <- do.call(
  rbind,
  lapply(panel_order, function(panel_tag) {
    data.frame(
      tag = factor(panel_tag, levels = panel_order),
      tag_label = panel_tag,
      stringsAsFactors = FALSE
    )
  })
)

figure_4 <- ggplot(facet_df, aes(x = comparison, y = feature_id)) +
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
    na.rm = TRUE,
    show.legend = FALSE
  ) +
  geom_text(
    data = tag_df,
    aes(x = -Inf, y = Inf, label = tag_label),
    inherit.aes = FALSE,
    hjust = 1.15,
    vjust = -0.4,
    size = 4.2,
    fontface = "bold"
  ) +
  facet_wrap(
    ~tag,
    ncol = 1,
    scales = "free_y",
    labeller = strip_labeller
  ) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL, caption = caption_text) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(colour = "black", size = 8, face = "bold"),
    legend.text = element_text(colour = "black", size = 7, face = "italic"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    plot.caption = element_text(face = "bold", size = 8),
    plot.margin = margin(t = 6, r = 4, b = 4, l = 24, unit = "pt")
  )

print(figure_4)

output_dir <- file.path(work_dir, "figures", "main_figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_file <- file.path(output_dir, "figure_4.pdf")

pdf(output_file, width = 14, height = 8, useDingbats = FALSE)
print(figure_4)
dev.off()
