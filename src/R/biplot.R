library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(reshape2)


# PERMANOVA results should be in the form "test_statistic,R_squared,p_value"
# Assuming permanova_results is passed as a parameter or loaded from elsewhere
# permanova_results should contain pseudo_F, r_squared, and p_value components
permanova_stats <- list(
  "test statistic" = as.numeric(permanova_results$pseudo_F),
  "R²" = as.numeric(permanova_results$r_squared),
  "p-value" = as.numeric(permanova_results$p_value)
)


# Pull out title
plot_title_names <- title_list[[2]]
# Prepare the title with stress value and PERMANOVA results
plot_title_stats <- paste0(
  "pseudo-F: ", round(permanova_stats[["test statistic"]], 6), ", ",
  "R²: ", round(permanova_stats[["R²"]], 6), ", ",
  "p-value: ", round(permanova_stats[["p-value"]], 4)
)


primary_variable <- title_list[[1]]



if (primary_variable == "Diagnosis") {
  ord_df <- ord_df %>%
    mutate(
      !!primary_variable := factor(
        .data[[primary_variable]],
        levels = c("Healthy", "GORD", "BO", "Dysplasia", "OAC", "Metastatic")
      )
    )
} else if (primary_variable == "sample_location") {
  ord_df <- ord_df %>%
    mutate(
      !!primary_variable := factor(
        .data[[primary_variable]],
        levels = c(
          "biopsy_location_1", "biopsy_location_2",
          "biopsy_location_3", "biopsy_location_4",
          "biopsy_location_5"
        )
      )
    )
} else {
  # Error handling for unexpected primary_variable values
  stop(paste("Unexpected primary_variable value:", primary_variable))
}




cb_palette <- c(
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#F0E442",
  "#0072B2",
  "#D55E00",
  "#CC79A7"
)


num_treatments <- length(levels(ord_df[[primary_variable]]))

print(num_treatments)

# Create the final color list with the first being grey
final_colors <- c("grey", cb_palette[1:(num_treatments - 1)])



# Convert primary_variable string to a formula dynamically
formula <- reformulate(primary_variable, response = "cbind(PC1, PC2)")

# Calculate centroids
centroids <- aggregate(formula, data = ord_df, mean)

# Merge data
rpca_data <- merge(ord_df, centroids,
  by = primary_variable,
  suffixes = c("", "_centroid")
)


# Use ggplot2 to create the biplot
biplot <- ggplot(data = rpca_data, aes(x = PC1, y = PC2, color = !!sym(primary_variable))) +
  geom_point(aes(fill = !!sym(primary_variable), shape = !!sym(primary_variable)),
    size = 3, stroke = 1.2
  ) +
  # Add ellipses by primary variable
  stat_ellipse(aes(fill = !!sym(primary_variable)),
    type = "norm",
    level = 0.95,
    alpha = 0.8,
    linetype = "dashed",
    linewidth = 1
  ) +
  # geom_point(aes(x = PC1_centroid, y = PC2_centroid, color = !!sym(primary_variable)), size = 4) +
  # geom_segment(aes(x = PC1_centroid, y = PC2_centroid, xend = PC1, yend = PC2, color = !!sym(primary_variable)), size = 1) +
  labs(
    title = paste0(plot_title_names, "\n", plot_title_stats),
    x = paste0("PC1 (", round(explained_variance[[1]] * 100, 2), "%)"),
    y = paste0("PC2 (", round(explained_variance[[2]] * 100, 2), "%)"),
  ) +
  geom_segment(
    data = top_feats_taxa_df, aes(x = 0, y = 0, xend = PC1 * 2, yend = PC2 * 2),
    arrow = arrow(length = unit(0.2, "cm")), color = "black"
  ) +
  geom_text_repel(
    data = top_feats_taxa_df, aes(x = PC1 * 2, y = PC2 * 2, label = label),
    color = "black", size = 3
  ) +
  scale_color_manual(values = final_colors) +
  scale_fill_manual(values = final_colors) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 8)) + # Add this line
  theme_classic() +
  theme(
    axis.line.y = element_line(linetype = 1, linewidth = 1.3, colour = "black"),
    axis.line.x = element_line(linetype = 1, linewidth = 1.3, colour = "black"),
    axis.title.y = element_text(face = "bold", size = 10),
    axis.title.x = element_text(face = "bold", size = 10),
    axis.text.x = element_text(face = "bold", size = 10, color = "black"),
    axis.text.y = element_text(face = "bold", size = 10, color = "black"),
    legend.position = "right",
    legend.title = element_text(colour = "black", size = 10, face = "bold"),
    legend.text = element_text(colour = "black", size = 8, face = "italic"),
    plot.tag = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10, hjust = 0.5, face = "italic", color = "black"),
    plot.title = element_text(size = 12, hjust = 0.5, face = "bold", color = "black")
  )
pdf(output_file, width = 14, height = 7.88, useDingbats = FALSE)
print(biplot)
dev.off()
