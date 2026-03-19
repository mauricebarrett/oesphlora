library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)


# Match working directory to analysis output location
work_dir <- "/home/mossy/projects/phd/oesphlora"

# Diagnosis of interest for Figure 5 (Barrett oesophagus)
diagnosis <- "BO"

# Path to alpha diversity comparison results between locations for this diagnosis
alpha_results_file <- file.path(
  work_dir,
  "diversity_metrics",
  "alpha_diversity",
  "results",
  "location",
  paste0("alpha_diversity_location_comparison_", diagnosis, ".csv")
)

alpha_df <- read.csv(alpha_results_file, stringsAsFactors = FALSE)

if (nrow(alpha_df) == 0) {
  stop("Alpha diversity results file is empty for diagnosis ", diagnosis, ".")
}

# Use z-score scaled differences for the heatmap
plot_column <- "zscore_diff"
legend_label <- "Z-score difference"
caption_text <- "* <= 0.05, ** <= 0.01, *** <= 0.001"

# Ensure required columns are present
required_cols <- c(
  "metric",
  "comparison",
  plot_column,
  "significance"
)

missing_cols <- setdiff(required_cols, colnames(alpha_df))
if (length(missing_cols) > 0) {
  stop(
    paste(
      "Missing required columns in alpha diversity results:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

# Tidy up comparison labels: drop 'biopsy_location_' prefix
alpha_df$comparison <- gsub("biopsy_location_", "", alpha_df$comparison)
alpha_df$comparison <- gsub("_", " ", alpha_df$comparison)

# Build matrix for symmetric colour scaling
data_matrix <- dcast(
  alpha_df,
  metric ~ comparison,
  value.var = plot_column
)
rownames(data_matrix) <- data_matrix$metric
data_matrix <- data_matrix[, -1, drop = FALSE]
data_matrix[is.na(data_matrix)] <- 0

# Set symmetric colour scale
max_abs_val <- max(abs(alpha_df[[plot_column]]), na.rm = TRUE)
if (!is.finite(max_abs_val) || max_abs_val == 0) {
  max_abs_val <- 1
}
color_limits <- c(-max_abs_val, max_abs_val)

# Order metrics as they appear in the CSV
alpha_df$metric <- factor(alpha_df$metric, levels = unique(alpha_df$metric))

# Order comparisons numerically by location index
alpha_df$comparison <- factor(
  alpha_df$comparison,
  levels = sort(unique(alpha_df$comparison))
)

figure_5 <- ggplot(
  alpha_df,
  aes(x = comparison, y = metric)
) +
  geom_tile(
    aes(fill = .data[[plot_column]]),
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "darkblue",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    space = "Lab",
    na.value = "grey90",
    name = legend_label,
    limits = color_limits
  ) +
  geom_text(
    aes(label = significance),
    color = "black",
    size = 4,
    vjust = 0.7,
    na.rm = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = NULL,
    caption = caption_text
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 12,
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 10,
      face = "bold",
      color = "black"
    ),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(
      colour = "black",
      size = 10,
      face = "bold"
    ),
    legend.text = element_text(
      colour = "black",
      size = 10,
      face = "italic"
    ),
    plot.caption = element_text(face = "bold")
  )

print(figure_5)

# Save main Figure 5
output_dir <- file.path(work_dir, "figures", "main_figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_file <- file.path(output_dir, "figure_5.pdf")

pdf(output_file, width = 14, height = 7.88, useDingbats = FALSE)
print(figure_5)
dev.off()

