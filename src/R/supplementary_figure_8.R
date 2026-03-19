library(ggplot2)


# Match working directory to analysis output location
work_dir <- "/home/mossy/projects/phd/oesphlora"

# Clinical categories to include in the supplementary alpha-diversity figure.
# BO is shown separately in main Figure 4.
diagnoses <- c("Healthy", "GORD", "Dysplasia", "OAC", "Metastatic")

caption_text <- "* <= 0.05, ** <= 0.01, *** <= 0.001"


get_results_file <- function(diagnosis) {
  file.path(
    work_dir,
    "diversity_metrics",
    "alpha_diversity",
    "results",
    "location",
    paste0("alpha_diversity_location_comparison_", diagnosis, ".csv")
  )
}


read_alpha_results <- function(diagnosis) {
  results_file <- get_results_file(diagnosis)

  if (!file.exists(results_file)) {
    stop(
      paste(
        "Alpha diversity location-comparison results file not found for diagnosis",
        diagnosis,
        "at:",
        results_file
      )
    )
  }

  alpha_df <- read.csv(results_file, stringsAsFactors = FALSE)

  if (nrow(alpha_df) == 0) {
    stop("Alpha diversity results file is empty for diagnosis ", diagnosis, ".")
  }

  required_cols <- c("metric", "comparison", "zscore_diff", "significance")
  missing_cols <- setdiff(required_cols, colnames(alpha_df))

  if (length(missing_cols) > 0) {
    stop(
      paste(
        "Missing required columns in alpha diversity results for",
        diagnosis,
        ":",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  alpha_df$comparison <- gsub("biopsy_location_", "", alpha_df$comparison)
  alpha_df$comparison <- gsub("_", " ", alpha_df$comparison)
  alpha_df$diagnosis <- diagnosis

  alpha_df
}


alpha_results_list <- lapply(diagnoses, read_alpha_results)
facet_df <- do.call(rbind, alpha_results_list)

comparison_levels <- sort(unique(facet_df$comparison))
metric_levels <- unique(facet_df$metric)

facet_df$comparison <- factor(facet_df$comparison, levels = comparison_levels)
facet_df$metric <- factor(facet_df$metric, levels = metric_levels)

panel_order <- LETTERS[seq_along(diagnoses)]
tag_lookup <- setNames(panel_order, diagnoses)
facet_df$tag <- factor(tag_lookup[facet_df$diagnosis], levels = panel_order)

strip_labeller <- as_labeller(setNames(diagnoses, panel_order))

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

max_abs_val <- max(abs(facet_df$zscore_diff), na.rm = TRUE)
if (!is.finite(max_abs_val) || max_abs_val == 0) {
  max_abs_val <- 1
}

supplementary_figure_8 <- ggplot(facet_df, aes(x = comparison, y = metric)) +
  geom_tile(aes(fill = zscore_diff), color = "black") +
  scale_fill_gradient2(
    low = "darkblue",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    space = "Lab",
    na.value = "grey90",
    name = "Z-score difference",
    limits = c(-max_abs_val, max_abs_val)
  ) +
  geom_text(
    aes(label = significance),
    color = "black",
    size = 3.3,
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
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(colour = "black", size = 10, face = "bold"),
    legend.text = element_text(colour = "black", size = 10, face = "italic"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    plot.caption = element_text(face = "bold", size = 8),
    plot.margin = margin(t = 6, r = 4, b = 4, l = 24, unit = "pt")
  )

print(supplementary_figure_8)

output_dir <- file.path(work_dir, "figures", "supplementary_figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_file <- file.path(output_dir, "supplementary_figure_8.pdf")

pdf(output_file, width = 8.27, height = 11.69, useDingbats = FALSE)
print(supplementary_figure_8)
dev.off()
