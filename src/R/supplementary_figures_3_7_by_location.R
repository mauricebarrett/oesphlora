library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)


# Set working directory (match analysis output location)
work_dir <- "/home/mossy/projects/phd/oesphlora"

# Base paths to DAA results (comparing diagnoses per biopsy location)
daa_diagnosis_dir <- file.path(
  work_dir,
  "differential_abundance_analysis",
  "diagnosis"
)

# Supplementary figures for all biopsy locations
locations <- paste0("biopsy_location_", c(1, 2, 3, 4, 5))

# Default ranks to include, in fixed panel order
default_ranks <- c("asv", "species", "genus")

# Output directory for supplementary figures
output_dir <- file.path(work_dir, "figures", "supplementary_figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# A4 page size in inches (portrait)
a4_width_in <- 8.27
a4_height_in <- 11.69

# Location to supplementary figure number mapping
figure_number_map <- c(
  biopsy_location_1 = 3,
  biopsy_location_2 = 4,
  biopsy_location_3 = 5,
  biopsy_location_4 = 6,
  biopsy_location_5 = 7
)


get_ranks_for_location <- function(location) {
  # For biopsy location 3, only include species and genus panels.
  if (location == "biopsy_location_3") {
    return(c("species", "genus"))
  }
  default_ranks
}


prepare_daa_heatmap_df <- function(daa_df) {
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

  if (nrow(data_matrix) > 1 && ncol(data_matrix) > 0) {
    hc <- hclust(dist(data_matrix, method = "euclidean"), method = "ward.D2")
    ordered_features <- rownames(data_matrix)[hc$order]
    plot_full_df$feature_id <- factor(plot_full_df$feature_id, levels = ordered_features)
  }

  plot_full_df$comparison_standard <- factor(
    plot_full_df$comparison,
    levels = all_comparisons
  )

  plot_full_df
}


resolve_results_file <- function(location, rank_label) {
  if (rank_label == "asv") {
    file.path(
      daa_diagnosis_dir,
      "asv_level",
      paste0(
        "differential_abundance_results_comparing_diagnoses_",
        location,
        "_filtered.csv"
      )
    )
  } else {
    file.path(
      daa_diagnosis_dir,
      paste0(rank_label, "_level"),
      paste0(
        "daa_",
        rank_label,
        "_results_comparing_diagnoses_",
        location,
        "_filtered.csv"
      )
    )
  }
}


# Returns prepared heatmap data for one panel, or NULL if unplottable.
# Caller adds panel_tag and rbinds for faceting.
get_panel_data <- function(location, rank_label) {
  results_file <- resolve_results_file(location, rank_label)

  if (!file.exists(results_file)) {
    message(
      "Omitting ",
      location,
      " ",
      rank_label,
      " panel: no file at ",
      results_file
    )
    return(NULL)
  }

  daa_df <- read.csv(results_file, stringsAsFactors = FALSE)
  required_cols <- c("feature_id", "comparison", "log2FoldChange", "padj")
  if (nrow(daa_df) == 0 || !all(required_cols %in% colnames(daa_df))) {
    message(
      "Omitting ",
      location,
      " ",
      rank_label,
      " panel: no plottable DAA rows"
    )
    return(NULL)
  }

  daa_df <- daa_df %>%
    filter(!is.na(feature_id), !is.na(comparison))

  if (nrow(daa_df) == 0 || length(unique(daa_df$feature_id)) < 2) {
    message(
      "Omitting ",
      location,
      " ",
      rank_label,
      " panel: insufficient features"
    )
    return(NULL)
  }

  if (rank_label == "asv") {
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
          if (is.na(val)) "unclassified" else val
        }
      )
      daa_df$feature_id <- paste0(daa_df$feature_id, "=", lowest_classified)
    }
  }

  plot_full_df <- prepare_daa_heatmap_df(daa_df)
  if (nrow(plot_full_df) == 0) {
    message(
      "Omitting ",
      location,
      " ",
      rank_label,
      " panel: empty plot data after prep"
    )
    return(NULL)
  }

  plot_full_df
}


make_location_figure <- function(location) {
  significance_caption <- "* <= 0.05, ** <= 0.01, *** <= 0.001"
  ranks <- get_ranks_for_location(location)
  panel_dfs <- list()
  panel_feature_levels <- character(0)

  for (rank_label in ranks) {
    panel_data <- get_panel_data(location = location, rank_label = rank_label)
    if (is.null(panel_data)) next
    panel_tag <- LETTERS[[length(panel_dfs) + 1]]
    panel_data$panel_tag <- panel_tag
    panel_data$rank_label <- rank_label
    panel_data$panel_feature <- paste(
      panel_tag,
      as.character(panel_data$feature_id),
      sep = "___"
    )
    panel_feature_levels <- c(
      panel_feature_levels,
      paste(panel_tag, levels(panel_data$feature_id), sep = "___")
    )
    panel_dfs[[length(panel_dfs) + 1]] <- panel_data
  }

  if (length(panel_dfs) == 0) {
    message("No plottable panels for ", location, "; skipping output.")
    return(invisible(NULL))
  }

  facet_df <- do.call(rbind, panel_dfs)
  facet_df$panel_feature <- factor(
    facet_df$panel_feature,
    levels = panel_feature_levels
  )
  panel_order <- LETTERS[seq_along(panel_dfs)]
  facet_df$panel_tag <- factor(facet_df$panel_tag, levels = panel_order)

  rank_labels <- c(
    asv = "ASV level",
    species = "Species level",
    genus = "Genus level"
  )
  strip_labeller <- ggplot2::as_labeller(
    setNames(
      rank_labels[ranks[seq_along(panel_dfs)]],
      panel_order
    )
  )

  tag_df <- data.frame(
    panel_tag = factor(panel_order, levels = panel_order),
    tag_label = panel_order,
    stringsAsFactors = FALSE
  )

  p <- ggplot(facet_df, aes(x = comparison_standard, y = panel_feature)) +
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
      size = 3.5,
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
      ~ panel_tag,
      ncol = 1,
      scales = "free_y",
      labeller = strip_labeller
    ) +
    scale_y_discrete(labels = function(x) sub("^[A-Z]___", "", x)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL, caption = significance_caption) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 6, face = "bold", color = "black"),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.title = element_text(colour = "black", size = 9, face = "bold"),
      legend.text = element_text(colour = "black", size = 8, face = "italic"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 9),
      plot.caption = element_text(face = "bold", size = 8),
      plot.margin = margin(t = 6, r = 4, b = 4, l = 24, unit = "pt")
    )

  figure_number <- unname(figure_number_map[location])
  if (is.na(figure_number)) {
    stop("No supplementary figure number mapped for ", location)
  }

  output_file <- file.path(
    output_dir,
    paste0("supplementary_figure_", figure_number, ".pdf")
  )

  pdf(output_file, width = a4_width_in, height = a4_height_in, useDingbats = FALSE)
  print(p)
  dev.off()

  message("Wrote ", output_file)
}


for (loc in locations) {
  make_location_figure(loc)
}
