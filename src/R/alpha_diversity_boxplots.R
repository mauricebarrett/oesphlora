library(ggplot2)
library(ggpubr)
library(tidyr)


# reorder levels of the primary variable
merged_df <- merged_df %>%
    mutate(
        !!primary_variable := factor(
            .data[[primary_variable]],
            levels = c("Healthy", "GORD", "BO", "Dysplasia", "OAC", "Metastatic")
        )
    )


# Create all possible pairwise comparisons without filtering
all_groups <- levels(factor(merged_df[[primary_variable]]))
pairwise_comparisons <- combn(all_groups, 2, simplify = FALSE)

# Print to verify comparisons
print("All pairwise comparisons:")
print(pairwise_comparisons)

# Define a colorblind-friendly palette
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
num_groups <- length(all_groups)
final_colors <- c("grey", cb_palette[1:min(num_groups-1, length(cb_palette))])

# Split the data by metric and create separate plots
metrics <- unique(merged_df$alpha_diversity_metric)
plot_list <- list()

for (metric in metrics) {
  metric_data <- merged_df[merged_df$alpha_diversity_metric == metric,]
  
  p <- ggplot(metric_data, aes(x = .data[[primary_variable]], y = value, fill = .data[[primary_variable]])) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values = final_colors) +
    ggtitle(metric) +
    theme(legend.position = "none") +
    stat_compare_means(
      comparisons = pairwise_comparisons,
      label = "p.signif", 
      method = "wilcox.test", 
      p.adjust.method = NULL,
      step.increase = 0.1
    ) +
    theme_classic() +
    labs(title = metric, x = NULL, y = NULL) +
    theme(
        legend.position = "none",
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        axis.line.y = element_line(linetype = 1, linewidth = 0.5, colour = "black"),
        axis.line.x = element_line(linetype = 1, linewidth = 0.5, colour = "black"),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 12, color = "black")
    )
  
  plot_list[[metric]] <- p
}

# Combine all plots
combined_plot <- ggpubr::ggarrange(plotlist = plot_list, ncol = 2, nrow = 2)

# Save the plot
temp_file <- tempfile(fileext = ".png")
png(temp_file, width = 16, height = 9, units = "in", res = 300)
print(combined_plot)
dev.off()

temp_file
