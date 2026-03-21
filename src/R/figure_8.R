library(ggplot2)
library(dplyr)
library(reshape2)


work_dir <- "/home/mossy/projects/phd/oesphlora"

picrust2_location_dir <- file.path(
  work_dir,
  "picrust2_differential_abundance_analysis",
  "location"
)

diagnosis <- "BO"

# Target max features per panel; actual count is min(this, available in each file)
# so every panel shows the same number of rows after filtering.
top_n_features <- 30L
# Rank features using either greatest |log2FC| first ("log2fc") or strongest
# significance first ("significance", i.e. smallest min padj across comparisons).
feature_rank_method <- "log2fc"


caption_text <- "* <= 0.05, ** <= 0.01, *** <= 0.001"
location_levels <- paste0("biopsy_location_", 1:5)
all_comparisons_raw <- combn(location_levels, 2, function(x) paste(x[1], "vs", x[2]))
all_comparisons <- gsub("biopsy_location_", "", all_comparisons_raw)


filter_daa_to_top_features <- function(daa_df, n, rank_method = "log2fc") {
  rank_method <- match.arg(rank_method, c("log2fc", "significance"))

  feat_stats <- daa_df %>%
    dplyr::group_by(feature_id) %>%
    dplyr::summarise(
      min_padj = min(padj, na.rm = TRUE),
      max_abs_lfc = max(abs(log2FoldChange), na.rm = TRUE),
      .groups = "drop"
    )

  feat_stats$min_padj[!is.finite(feat_stats$min_padj)] <- 1
  feat_stats$max_abs_lfc[!is.finite(feat_stats$max_abs_lfc)] <- 0

  if (rank_method == "log2fc") {
    feat_stats <- feat_stats %>%
      dplyr::arrange(dplyr::desc(max_abs_lfc), min_padj)
  } else {
    feat_stats <- feat_stats %>%
      dplyr::arrange(min_padj, dplyr::desc(max_abs_lfc))
  }

  n_keep <- min(as.integer(n), nrow(feat_stats))
  keep_ids <- head(feat_stats$feature_id, n_keep)

  daa_df %>% dplyr::filter(feature_id %in% keep_ids)
}


read_daa_csv <- function(pathway_name, subfolder, file_prefix) {
  daa_file <- file.path(
    picrust2_location_dir,
    subfolder,
    paste0(
      "daa_",
      file_prefix,
      "_results_comparing_locations_",
      diagnosis,
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

  daa_df %>%
    dplyr::filter(!is.na(feature_id), !is.na(comparison))
}


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
  data_matrix <- data_matrix[, -1, drop = FALSE]
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


# 1) Load long tables (no expand.grid / heatmap prep yet).
ko_raw <- read_daa_csv("KEGG KO", "kegg_level", "kegg")
ec_raw <- read_daa_csv("Enzyme Commission (EC)", "ec_level", "ec")
metacyc_raw <- read_daa_csv("MetaCyc", "metacyc_level", "metacyc")

n_ko <- dplyr::n_distinct(ko_raw$feature_id)
n_ec <- dplyr::n_distinct(ec_raw$feature_id)
n_mc <- dplyr::n_distinct(metacyc_raw$feature_id)

if (min(n_ko, n_ec, n_mc) < 1L) {
  stop("At least one pathway file has no features after filtering NA rows.")
}

# Same feature count in every panel: cap by the smallest ontology and top_n_features.
n_per_panel <- min(as.integer(top_n_features), n_ko, n_ec, n_mc)

# 2) Filter each long table to exactly n_per_panel features (ranked within ontology).
ko_f <- filter_daa_to_top_features(
  ko_raw,
  n = n_per_panel,
  rank_method = feature_rank_method
)
ec_f <- filter_daa_to_top_features(
  ec_raw,
  n = n_per_panel,
  rank_method = feature_rank_method
)
metacyc_f <- filter_daa_to_top_features(
  metacyc_raw,
  n = n_per_panel,
  rank_method = feature_rank_method
)

stopifnot(
  dplyr::n_distinct(ko_f$feature_id) == n_per_panel,
  dplyr::n_distinct(ec_f$feature_id) == n_per_panel,
  dplyr::n_distinct(metacyc_f$feature_id) == n_per_panel
)

# 3) Build heatmap-ready rows, then combine panels.
ko_df <- prepare_daa_heatmap_df(ko_f)
ec_df <- prepare_daa_heatmap_df(ec_f)
metacyc_df <- prepare_daa_heatmap_df(metacyc_f)

facet_df <- rbind(
  cbind(ko_df, tag = "A"),
  cbind(ec_df, tag = "B"),
  cbind(metacyc_df, tag = "C")
)

panel_order <- c("A", "B", "C")
facet_df$tag <- factor(facet_df$tag, levels = panel_order)
strip_labeller <- as_labeller(
  c(
    A = "KEGG Orthologues",
    B = "Enzyme Commission",
    C = "MetaCyc pathways"
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

figure_8 <- ggplot(facet_df, aes(x = comparison, y = feature_id)) +
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
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 6, face = "bold", color = "black"),
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

print(figure_8)

output_dir <- file.path(work_dir, "figures", "main_figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_file <- file.path(output_dir, "figure_8.pdf")



pdf(
  output_file,
  width = 14,
  height = 8,
  useDingbats = FALSE
)
print(figure_8)
dev.off()
