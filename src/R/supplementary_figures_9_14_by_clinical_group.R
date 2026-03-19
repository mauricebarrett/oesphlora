library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)


# Set working directory (match analysis output location)
work_dir <- "/home/mossy/projects/phd/oesphlora"

# Base paths to DAA results (comparing biopsy locations per clinical group)
daa_location_dir <- file.path(
  work_dir,
  "differential_abundance_analysis",
  "location"
)

# Clinical groups to include in supplementary figures
clinical_groups <- c("Healthy", "GORD", "BO", "Dysplasia", "OAC", "Metastatic")

# Default ranks to include, in fixed panel order
default_ranks <- c("asv", "species", "genus")

# Output directory for supplementary figures
output_dir <- file.path(work_dir, "figures", "supplementary_figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# A4 page size in inches (portrait)
a4_width_in <- 8.27
a4_height_in <- 11.69

# Clinical group to supplementary figure number mapping
figure_number_map <- c(
  Healthy = 9,
  GORD = 10,
  BO = 11,
  Dysplasia = 12,
  OAC = 13,
  Metastatic = 14
)


get_ranks_for_clinical_group <- function(clinical_group) {
  # BO ASV-level heatmap is already shown in main Figure 7.
  if (clinical_group == "BO") {
    return(c("species", "genus"))
  }

  default_ranks
}


prepare_daa_heatmap_df <- function(daa_df) {
  daa_df$comparison <- gsub("biopsy_location_", "", daa_df$comparison)
  location_levels <- as.character(1:5)
  all_comparisons <- combn(location_levels, 2, function(x) paste(x[1], "vs", x[2]))

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
    dplyr::mutate(
      star = dplyr::case_when(
        padj < 0.001 ~ "***",
        padj < 0.01 ~ "**",
        padj < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::select(feature_id, comparison, log2FoldChange, star, padj)

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


resolve_results_file <- function(clinical_group, rank_label) {
  if (rank_label == "asv") {
    filtered_path <- file.path(
      daa_location_dir,
      "asv_level",
      paste0("differential_abundance_results_location_", clinical_group, "_filtered.csv")
    )
    unfiltered_path <- file.path(
      daa_location_dir,
      "asv_level",
      paste0("differential_abundance_results_location_", clinical_group, ".csv")
    )
  } else {
    filtered_path <- file.path(
      daa_location_dir,
      paste0(rank_label, "_level"),
      paste0("daa_", rank_label, "_results_location_", clinical_group, "_filtered.csv")
    )
    unfiltered_path <- file.path(
      daa_location_dir,
      paste0(rank_label, "_level"),
      paste0("daa_", rank_label, "_results_location_", clinical_group, ".csv")
    )
  }

  if (file.exists(filtered_path)) {
    return(filtered_path)
  }

  if (file.exists(unfiltered_path)) {
    return(unfiltered_path)
  }

  NA_character_
}


# Returns prepared heatmap data for one panel, or NULL if unplottable.
# Caller adds panel_tag and rbinds for faceting.
get_panel_data <- function(clinical_group, rank_label) {
  results_file <- resolve_results_file(clinical_group, rank_label)

  if (is.na(results_file)) {
    message(
      "Omitting ",
      clinical_group,
      " ",
      rank_label,
      " panel: no results file found"
    )
    return(NULL)
  }

  daa_df <- read.csv(results_file, stringsAsFactors = FALSE)
  required_cols <- c("feature_id", "comparison", "log2FoldChange", "padj")
  if (nrow(daa_df) == 0 || !all(required_cols %in% colnames(daa_df))) {
    message(
      "Omitting ",
      clinical_group,
      " ",
      rank_label,
      " panel: no plottable DAA rows"
    )
    return(NULL)
  }

  # If this is an unfiltered export, retain significant rows to keep the heatmap readable.
  if (!grepl("_filtered\\.csv$", results_file)) {
    daa_df <- daa_df %>% dplyr::filter(!is.na(padj), padj < 0.05)
  }

  daa_df <- daa_df %>%
    dplyr::filter(!is.na(feature_id), !is.na(comparison))

  if (nrow(daa_df) == 0 || length(unique(daa_df$feature_id)) < 2) {
    message(
      "Omitting ",
      clinical_group,
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
      clinical_group,
      " ",
      rank_label,
      " panel: empty plot data after prep"
    )
    return(NULL)
  }

  plot_full_df
}


make_clinical_group_figure <- function(clinical_group) {
  significance_caption <- "* <= 0.05, ** <= 0.01, *** <= 0.001"
  ranks <- get_ranks_for_clinical_group(clinical_group)
  panel_dfs <- list()
  panel_feature_levels <- character(0)

  for (rank_label in ranks) {
    panel_data <- get_panel_data(
      clinical_group = clinical_group,
      rank_label = rank_label
    )
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
    message("No plottable panels for ", clinical_group, "; skipping output.")
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
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6, face = "bold"),
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

  figure_number <- unname(figure_number_map[clinical_group])
  if (is.na(figure_number)) {
    stop("No supplementary figure number mapped for ", clinical_group)
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


for (clinical_group in clinical_groups) {
  make_clinical_group_figure(clinical_group)
}
