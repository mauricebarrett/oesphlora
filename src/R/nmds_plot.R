#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)



# PERMANOVA results should be in the form "test_statistic,R_squared,p_value"
permanova_results <- list(
  "test statistic" = as.numeric(permanova_results$pseudo_F),
  "R²" = as.numeric(permanova_results$r_squared),
  "p-value" = as.numeric(permanova_results$p_value)
)




# Pull out title
plot_title_names <- title_list[[2]]


# Prepare the title with stress value and PERMANOVA results
plot_title_stats <- paste0(
  "Stress: ", round(stress_value, 6), ", ",
  "PERMANOVA pseudo-F: ", round(permanova_results[["test statistic"]], 6), ", ",
  "R²: ", round(permanova_results[["R²"]], 6), ", ",
  "p-value: ", round(permanova_results[["p-value"]], 4)
)




# Add row names as a column manually for both data frames
nmds_coordinates_df$SampleID <- rownames(nmds_coordinates_df)
metadata$SampleID <- rownames(metadata)

# Perform the join based on 'SampleID'
nmds_data <- merge(nmds_coordinates_df, metadata, by = "SampleID")


primary_variable <- title_list[[1]]


if (primary_variable == "Diagnosis") {
  nmds_data <- nmds_data %>%
    mutate(
      !!primary_variable := factor(
        .data[[primary_variable]],
        levels = c("Healthy", "GORD", "BO", "Dysplasia", "OAC", "Metastatic")
      )
    )
} else if (primary_variable == "sample_location") {
  nmds_data <- nmds_data %>%
    mutate(
      !!primary_variable := factor(
        .data[[primary_variable]],
        levels = c("biopsy_location_1", "biopsy_location_2", "biopsy_location_3", "biopsy_location_4", "biopsy_location_5")
      )
    )
} else {
  # Error handling for unexpected primary_variable values
  stop(paste("Unexpected primary_variable value:", primary_variable))
}




# Define a colorblind-friendly palette (Okabe-Ito or similar)
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Ensure the color vector length matches `treatments_used`
num_treatments <- length(levels(nmds_data[[primary_variable]]))

print(num_treatments)

# Create the final color list with the first being grey
final_colors <- c("grey", cb_palette[1:(num_treatments - 1)])


# Convert primary_variable string to a formula dynamically
formula <- reformulate(primary_variable, response = "cbind(NMDS1, NMDS2)")

# Calculate centroids
centroids <- aggregate(formula, data = nmds_data, mean)

# Merge data
bc_data <- merge(nmds_data,
  centroids,
  by = primary_variable,
  suffixes = c("", "_centroid")
)


# Create the NMDS plot
nmds_plot <- ggplot(bc_data, aes(x = NMDS1, y = NMDS2, color = !!sym(primary_variable))) +
  geom_point(aes(fill = !!sym(primary_variable), shape = !!sym(primary_variable)),
    size = 3, stroke = 1.2
  ) +
  geom_point(aes(x = NMDS1_centroid, y = NMDS2_centroid, color = !!sym(primary_variable)), size = 4) +
  stat_ellipse(aes(fill = !!sym(primary_variable)),
    type = "norm",
    level = 0.95,
    alpha = 0.8,
    linetype = "dashed",
    linewidth = 1
  ) +
  # stat_ellipse(geom = "polygon", aes(x = NMDS1, y = NMDS2, group = primary_variable, color = primary_variable, fill = primary_variable), level = 0.95, type = "norm", size = 1, alpha = 0.3) +
  labs(
    title = paste0(plot_title_names, "\n", plot_title_stats),
    x = "NMDS1",
    y = "NMDS2"
  ) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 8)) + # Add this line
  scale_color_manual(values = final_colors) +
  scale_fill_manual(values = final_colors) +
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

# Save the plot to a file
pdf(output_file, width = 14, height = 8, useDingbats = FALSE)
print(nmds_plot)
dev.off()
