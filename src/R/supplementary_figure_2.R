library(ggplot2)


# Match working directory to analysis output location
work_dir <- "/home/mossy/projects/phd/oesphlora"

caption_text <- "* <= 0.05, ** <= 0.01, *** <= 0.001"
group_levels <- c("Healthy", "GORD", "BO", "Dysplasia", "OAC", "Metastatic")
all_comparisons <- combn(group_levels, 2, function(x) paste(x[1], "vs", x[2]))


alpha_results_dir <- file.path(
  work_dir,
  "diversity_metrics",
  "alpha_diversity",
  "results",
  "diagnosis"
)

# Read each biopsy-site results table explicitly (one-by-one)
alpha_location_1 <- read.csv(
  file.path(
    alpha_results_dir,
    "alpha_diversity_diagnosis_comparison_biopsy_location_1.csv"
  ),
  stringsAsFactors = FALSE
)

alpha_location_2 <- read.csv(
  file.path(
    alpha_results_dir,
    "alpha_diversity_diagnosis_comparison_biopsy_location_2.csv"
  ),
  stringsAsFactors = FALSE
)

alpha_location_3 <- read.csv(
  file.path(
    alpha_results_dir,
    "alpha_diversity_diagnosis_comparison_biopsy_location_3.csv"
  ),
  stringsAsFactors = FALSE
)

alpha_location_4 <- read.csv(
  file.path(
    alpha_results_dir,
    "alpha_diversity_diagnosis_comparison_biopsy_location_4.csv"
  ),
  stringsAsFactors = FALSE
)

alpha_location_5 <- read.csv(
  file.path(
    alpha_results_dir,
    "alpha_diversity_diagnosis_comparison_biopsy_location_5.csv"
  ),
  stringsAsFactors = FALSE
)


# rbind 
facet_df <- rbind(
  cbind(alpha_location_1, location_label = "Biopsy location 1", tag = "A"),
  cbind(alpha_location_2, location_label = "Biopsy location 2", tag = "B"),
  cbind(alpha_location_3, location_label = "Biopsy location 3", tag = "C"),
  cbind(alpha_location_4, location_label = "Biopsy location 4", tag = "D"),
  cbind(alpha_location_5, location_label = "Biopsy location 5", tag = "E")
)

facet_df$comparison <- gsub("_vs_", " vs ", facet_df$comparison)
facet_df$comparison <- factor(facet_df$comparison, levels = all_comparisons)

panel_order <- c("A", "B", "C", "D", "E")
facet_df$tag <- factor(facet_df$tag, levels = panel_order)
strip_labeller <- as_labeller(
  c(
    A = "Biopsy location 1",
    B = "Biopsy location 2",
    C = "Biopsy location 3",
    D = "Biopsy location 4",
    E = "Biopsy location 5"
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

supplementary_figure_2 <- ggplot(facet_df, aes(x = comparison, y = metric)) +
  geom_tile(aes(fill = zscore_diff), color = "black") +
  scale_fill_gradient2(
    low = "darkblue",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    space = "Lab",
    na.value = "grey90",
    name = "Z-score difference"
  ) +
  geom_text(
    aes(label = significance),
    color = "black",
    size = 5,
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
    legend.text = element_text(colour = "black", size = 8, face = "italic"),
    # legend.key.size = unit(8, "pt"),
    # legend.spacing.y = unit(2, "pt"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    plot.caption = element_text(face = "bold", size = 8),
    plot.margin = margin(t = 6, r = 4, b = 4, l = 24, unit = "pt")
  )

print(supplementary_figure_2)

output_dir <- file.path(work_dir, "figures", "supplementary_figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

output_file <- file.path(output_dir, "supplementary_figure_2.pdf")

pdf(output_file, width = 8.27, height = 11.69, useDingbats = FALSE)
print(supplementary_figure_2)
dev.off()
