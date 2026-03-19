library(ggplot2)
library(dplyr)
library(tidyr)
library(MASS)
library(cowplot)


# Match working directory to analysis output location
work_dir <- "/home/mossy/projects/phd/oesphlora"

bray_curtis_dir <- file.path(
  work_dir,
  "diversity_metrics",
  "beta_diversity",
  "location",
  "bray_curtis"
)

metadata_file <- file.path(work_dir, "metadata", "merged_data.csv")

if (!file.exists(metadata_file)) {
  stop(paste("Metadata file not found at:", metadata_file))
}

metadata_df <- read.csv(metadata_file, row.names = 1, stringsAsFactors = FALSE)
metadata_df$SampleID <- rownames(metadata_df)

diagnosis_levels <- c("Healthy", "GORD", "BO", "Dysplasia", "OAC", "Metastatic")
location_levels <- paste0("biopsy_location_", 1:5)
have_vegan <- requireNamespace("vegan", quietly = TRUE)

load_distance_matrix <- function(diagnosis) {
  distance_file <- file.path(
    bray_curtis_dir,
    paste0("bray_curtis_distance_matrix_comparing_locations_", diagnosis, ".csv")
  )

  if (!file.exists(distance_file)) {
    return(NULL)
  }

  distance_df <- read.csv(
    distance_file,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  distance_matrix <- as.matrix(distance_df)
  storage.mode(distance_matrix) <- "numeric"
  distance_matrix <- (distance_matrix + t(distance_matrix)) / 2
  diag(distance_matrix) <- 0
  distance_matrix
}

load_nmds_data <- function(diagnosis) {
  nmds_file <- file.path(
    bray_curtis_dir,
    paste0("bray_curtis_nmds_coordinates_", diagnosis, ".csv")
  )

  if (!file.exists(nmds_file)) {
    distance_matrix <- load_distance_matrix(diagnosis)

    if (is.null(distance_matrix)) {
      return(NULL)
    }

    if (nrow(distance_matrix) < 2) {
      return(NULL)
    }

    start_coords <- cmdscale(as.dist(distance_matrix), k = 2)

    nmds_fit <- tryCatch(
      isoMDS(as.dist(distance_matrix), y = start_coords, k = 2, trace = FALSE),
      error = function(e) NULL
    )

    if (is.null(nmds_fit)) {
      coords <- start_coords
    } else {
      coords <- nmds_fit$points
    }

    coords_df <- as.data.frame(coords, stringsAsFactors = FALSE)
    colnames(coords_df) <- c("NMDS1", "NMDS2")
    write.csv(coords_df, nmds_file, quote = FALSE)

    coords_df$SampleID <- rownames(coords_df)
    coords_df$Diagnosis <- diagnosis
    return(coords_df)
  }

  nmds_df <- read.csv(nmds_file, row.names = 1, stringsAsFactors = FALSE)
  nmds_df$SampleID <- rownames(nmds_df)
  nmds_df$Diagnosis <- diagnosis
  nmds_df
}

load_permanova_data <- function(diagnosis) {
  permanova_file <- file.path(
    bray_curtis_dir,
    paste0("bray_curtis_permanova_results_", diagnosis, ".csv")
  )

  if (!file.exists(permanova_file)) {
    if (!have_vegan) {
      return(NULL)
    }

    distance_matrix <- load_distance_matrix(diagnosis)

    if (is.null(distance_matrix)) {
      return(NULL)
    }

    sample_ids <- rownames(distance_matrix)
    metadata_subset <- metadata_df[sample_ids, , drop = FALSE]

    metadata_subset <- metadata_subset %>%
      filter(
        !is.na(sample_location),
        !is.na(patient_id)
      )

    if (nrow(metadata_subset) < 2 || dplyr::n_distinct(metadata_subset$sample_location) < 2) {
      return(NULL)
    }

    distance_matrix <- distance_matrix[
      metadata_subset$SampleID,
      metadata_subset$SampleID,
      drop = FALSE
    ]

    distance_dist <- as.dist(distance_matrix)
    adonis_result <- vegan::adonis2(
      distance_dist ~ sample_location,
      data = metadata_subset,
      permutations = 999,
      strata = metadata_subset[["patient_id"]]
    )

    pseudo_F <- adonis_result$F[1]
    p_value <- adonis_result$`Pr(>F)`[1]
    r_squared <- adonis_result$R2[1]

    dispersion_result <- vegan::betadisper(distance_dist, metadata_subset$sample_location)
    dispersion_test <- vegan::permutest(dispersion_result, permutations = 999)

    permanova_df <- data.frame(
      metric = c(
        "pseudo_F",
        "p_value",
        "r_squared",
        "dispersion_F",
        "dispersion_p_value"
      ),
      value = c(
        pseudo_F,
        p_value,
        r_squared,
        dispersion_test$tab$F[1],
        dispersion_test$tab$`Pr(>F)`[1]
      ),
      stringsAsFactors = FALSE
    )

    write.csv(permanova_df, permanova_file, row.names = FALSE, quote = FALSE)
    permanova_df$Diagnosis <- diagnosis
    return(permanova_df)
  }

  permanova_df <- read.csv(permanova_file, stringsAsFactors = FALSE)

  if (!all(c("metric", "value") %in% colnames(permanova_df))) {
    return(NULL)
  }

  permanova_df$Diagnosis <- diagnosis
  permanova_df
}

nmds_list <- lapply(diagnosis_levels, load_nmds_data)
nmds_list <- nmds_list[!vapply(nmds_list, is.null, logical(1))]

if (length(nmds_list) == 0) {
  stop(
    paste(
      "No Bray-Curtis NMDS coordinate files found in",
      bray_curtis_dir
    )
  )
}

permanova_list <- lapply(diagnosis_levels, load_permanova_data)
permanova_list <- permanova_list[!vapply(permanova_list, is.null, logical(1))]

nmds_data <- bind_rows(nmds_list)

plot_df <- nmds_data %>%
  inner_join(
    metadata_df %>% dplyr::select(SampleID, Diagnosis, sample_location),
    by = c("SampleID", "Diagnosis")
  )

if (nrow(plot_df) == 0) {
  stop("No samples matched between NMDS coordinates and metadata.")
}

plot_df$Diagnosis <- factor(plot_df$Diagnosis, levels = diagnosis_levels)
plot_df$sample_location <- factor(
  plot_df$sample_location,
  levels = location_levels
)

location_labels <- c(
  "Biopsy location 1",
  "Biopsy location 2",
  "Biopsy location 3",
  "Biopsy location 4",
  "Biopsy location 5"
)

plot_df$sample_location_label <- factor(
  c(
    "Biopsy location 1",
    "Biopsy location 2",
    "Biopsy location 3",
    "Biopsy location 4",
    "Biopsy location 5"
  )[as.integer(plot_df$sample_location)],
  levels = location_labels
)

permanova_labels <- NULL

if (length(permanova_list) > 0) {
  permanova_labels <- bind_rows(permanova_list) %>%
    mutate(value = as.numeric(value)) %>%
    pivot_wider(names_from = metric, values_from = value) %>%
    filter(Diagnosis %in% unique(as.character(plot_df$Diagnosis))) %>%
    mutate(
      Diagnosis = factor(Diagnosis, levels = diagnosis_levels),
      label = paste0(
        "PERMANOVA: pseudo-F: ", round(pseudo_F, 3),
        ", R²: ", round(r_squared, 3),
        ", p: ", round(p_value, 3)
      )
    )
}

panel_tags <- c("A", "B", "C", "D", "E", "F")
panel_label_df <- data.frame(
  Diagnosis = diagnosis_levels,
  panel_tag = panel_tags,
  stringsAsFactors = FALSE
)

if (!is.null(permanova_labels)) {
  panel_label_df <- panel_label_df %>%
    dplyr::left_join(
      permanova_labels %>% dplyr::select(Diagnosis, label),
      by = "Diagnosis"
    )
} else {
  panel_label_df$label <- NA_character_
}

location_palette <- c(
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#0072B2",
  "#D55E00"
)

build_diagnosis_plot <- function(diagnosis_name, panel_tag, permanova_label, show_legend = FALSE) {
  diagnosis_data <- plot_df %>%
    dplyr::filter(Diagnosis == diagnosis_name)
  subtitle_text <- if (is.na(permanova_label)) "" else permanova_label

  ggplot(
    diagnosis_data,
    aes(x = NMDS1, y = NMDS2, color = sample_location_label)
  ) +
    geom_point(
      aes(fill = sample_location_label, shape = sample_location_label),
      size = 1.8,
      stroke = 1.1
    ) +
    stat_ellipse(
      aes(group = sample_location_label),
      type = "norm",
      level = 0.5,
      alpha = 0.8,
      linetype = "dashed",
      linewidth = 1
    ) +
    scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
    scale_color_manual(values = location_palette, drop = FALSE) +
    scale_fill_manual(values = location_palette, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
    labs(
      x = "NMDS1",
      y = "NMDS2",
      color = "Biopsy location",
      fill = "Biopsy location",
      shape = "Biopsy location",
      tag = panel_tag,
      title = diagnosis_name,
      subtitle = subtitle_text
    ) +
    theme_classic() +
    theme(
      axis.line.y = element_line(linetype = 1, linewidth = 1.3, colour = "black"),
      axis.line.x = element_line(linetype = 1, linewidth = 1.3, colour = "black"),
      axis.title.y = element_text(face = "bold", size = 12),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.text.x = element_text(face = "bold", size = 10, color = "black"),
      axis.text.y = element_text(face = "bold", size = 10, color = "black"),
      legend.position = if (show_legend) "bottom" else "none",
      legend.background = element_rect(fill = "white", color = NA),
      legend.title = element_text(colour = "black", size = 11, face = "bold"),
      legend.text = element_text(colour = "black", size = 9, face = "italic"),
      plot.tag = element_text(face = "bold", size = 14),
      plot.tag.position = c(0.02, 0.98),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(face = "italic", size = 9, hjust = 0.5),
      plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
    )
}

plot_list <- lapply(seq_len(nrow(panel_label_df)), function(i) {
  build_diagnosis_plot(
    diagnosis_name = panel_label_df$Diagnosis[i],
    panel_tag = panel_label_df$panel_tag[i],
    permanova_label = panel_label_df$label[i],
    show_legend = FALSE
  )
})

legend_plot <- build_diagnosis_plot(
  diagnosis_name = panel_label_df$Diagnosis[1],
  panel_tag = panel_label_df$panel_tag[1],
  permanova_label = panel_label_df$label[1],
  show_legend = TRUE
)

shared_legend <- cowplot::get_legend(legend_plot)
panel_grid <- cowplot::plot_grid(plotlist = plot_list, ncol = 2, align = "hv")
figure_5 <- cowplot::plot_grid(panel_grid, shared_legend, ncol = 1, rel_heights = c(1, 0.12))

print(figure_5)

output_dir <- file.path(work_dir, "figures", "main_figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_file <- file.path(output_dir, "figure_6.pdf")

pdf(output_file, width = 14, height = 10, useDingbats = FALSE)
print(figure_5)
dev.off()
