library(ggplot2)
library(ggrepel)

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
    "pseudo-F: ", round(permanova_results[["test statistic"]], 6), ", ",
    "R²: ", round(permanova_results[["R²"]], 6), ", ",
    "p-value: ", round(permanova_results[["p-value"]], 4)
)


primary_variable <- title_list[[1]]


ord_df <- ord_df %>%
    mutate(
        !!primary_variable := factor(
            .data[[primary_variable]],
            levels = c("Healthy", "GORD", "BO", "Dysplasia", "OAC", "Metastatic")
        )
    )


# Define a colorblind-friendly palette (Okabe-Ito or similar)
cb_palette <- c(
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"
)

# Ensure the color vector length matches the number of levels in the primary variable
num_treatments <- length(levels(ord_df[[primary_variable]]))

print(num_treatments)

# Create the final color list with the first being grey
final_colors <- c("grey", cb_palette[1:(num_treatments - 1)])



# Convert primary_variable string to a formula dynamically
formula <- reformulate(primary_variable, response = "cbind(PC1, PC2)")

# Calculate centroids
centroids <- aggregate(formula, data = ord_df, mean)

# Merge data
rpca_data <- merge(ord_df,
    centroids,
    by = primary_variable,
    suffixes = c("", "_centroid")
)


biplot <- ggplot(data = rpca_data, aes(x = PC1, y = PC2, color = !!sym(primary_variable))) +
    geom_point(size = 2) +
    # geom_point(aes(x = PC1_centroid, y = PC2_centroid, color = !!sym(primary_variable)), size = 4) +
    # geom_segment(aes(x = PC1_centroid, y = PC2_centroid, xend = PC1, yend = PC2, color = !!sym(primary_variable)), size = 1) +
    labs(
        title = paste0(plot_title_names, "\n", plot_title_stats),
        x = paste0("PC1 (", round(explained_variance[[1]] * 100, 2), "%)"),
        y = paste0("PC2 (", round(explained_variance[[2]] * 100, 2), "%)"),
    ) +
    geom_segment(
        data = top_feats_taxa_df, aes(x = 0, y = 0, xend = PC1 * 1, yend = PC2 * 1),
        arrow = arrow(length = unit(0.2, "cm")), color = "black"
    ) +
    geom_text_repel(
        data = top_feats_taxa_df, aes(x = PC1 * 1, y = PC2 * 1, label = label),
        color = "black", size = 3
    ) +
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

temp_file <- tempfile(fileext = ".png")

png(temp_file, width = 16, height = 9, units = "in", res = 300)
print(biplot)
dev.off()

temp_file
