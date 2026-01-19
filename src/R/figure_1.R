library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
suppressPackageStartupMessages(library(ggthemes))


work_dir <- "/home/maurice/projects/phd/oesphlora"
bray_curtis_dir <- file.path(work_dir, "diversity_metrics", "beta_diversity", "diagnosis", "bray_curtis")

# read in the data
nmds_data_location_1 <- as.data.frame(read.csv(file.path(bray_curtis_dir, "bray_curtis_nmds_coordinates_biopsy_location_1.csv"), row.names = 1))
nmds_data_location_1$location <- "biopsy_location_1"
nmds_data_location_2 <- as.data.frame(read.csv(file.path(bray_curtis_dir, "bray_curtis_nmds_coordinates_biopsy_location_2.csv"), row.names = 1))
nmds_data_location_2$location <- "biopsy_location_2"
nmds_data_location_3 <- as.data.frame(read.csv(file.path(bray_curtis_dir, "bray_curtis_nmds_coordinates_biopsy_location_3.csv"), row.names = 1))
nmds_data_location_3$location <- "biopsy_location_3"
nmds_data_location_4 <- as.data.frame(read.csv(file.path(bray_curtis_dir, "bray_curtis_nmds_coordinates_biopsy_location_4.csv"), row.names = 1))
nmds_data_location_4$location <- "biopsy_location_4"
nmds_data_location_5 <- as.data.frame(read.csv(file.path(bray_curtis_dir, "bray_curtis_nmds_coordinates_biopsy_location_5.csv"), row.names = 1))
nmds_data_location_5$location <- "biopsy_location_5"


nmds_data <- rbind(nmds_data_location_1, nmds_data_location_2, nmds_data_location_3, nmds_data_location_4, nmds_data_location_5)


# # read in the permanova results
permanova_results_location_1 <- read.csv(file.path(bray_curtis_dir, "bray_curtis_permanova_results_biopsy_location_1.csv"))
permanova_results_location_1$location <- "biopsy_location_1"
permanova_results_location_2 <- read.csv(file.path(bray_curtis_dir, "bray_curtis_permanova_results_biopsy_location_2.csv"))
permanova_results_location_2$location <- "biopsy_location_2"
permanova_results_location_3 <- read.csv(file.path(bray_curtis_dir, "bray_curtis_permanova_results_biopsy_location_3.csv"))
permanova_results_location_3$location <- "biopsy_location_3"
permanova_results_location_4 <- read.csv(file.path(bray_curtis_dir, "bray_curtis_permanova_results_biopsy_location_4.csv"))
permanova_results_location_4$location <- "biopsy_location_4"
permanova_results_location_5 <- read.csv(file.path(bray_curtis_dir, "bray_curtis_permanova_results_biopsy_location_5.csv"))
permanova_results_location_5$location <- "biopsy_location_5"

# combine the permanova results
permanova_results <- rbind(permanova_results_location_1, permanova_results_location_2, permanova_results_location_3, permanova_results_location_4, permanova_results_location_5)

# Pivot wider to get pseudo_F, r_squared, and p_value as columns
permanova_results_wide <- pivot_wider(permanova_results, names_from = metric, values_from = value)

# Clean up location names in permanova results to match plot data
permanova_results_wide$location <- gsub("_", " ", permanova_results_wide$location)
permanova_results_wide$location <- gsub("biopsy", "Biopsy", permanova_results_wide$location)

# Create label with stats
permanova_results_wide$label <- paste0(
  "PERMANOVA: ",
  "pseudo-F: ", round(permanova_results_wide$pseudo_F, 3), ", ",
  "R²: ", round(permanova_results_wide$r_squared, 3), ", ",
  "p: ", round(permanova_results_wide$p_value, 3)
)

metadata_df <- read.csv(file.path("/home/maurice/projects/phd/oesphlora/metadata/merged_data.csv"), row.names = 1)

# merge the nmds data with the metadata
nmds_plot_data <- merge(nmds_data, metadata_df, by = "row.names")


nmds_plot_data$Diagnosis <- factor(nmds_plot_data$Diagnosis, levels = c("Healthy", "GORD", "BO", "Dysplasia", "OAC", "Metastatic"))
nmds_plot_data$location <- factor(nmds_plot_data$location, levels = c("biopsy_location_1", "biopsy_location_2", "biopsy_location_3", "biopsy_location_4", "biopsy_location_5"))
nmds_plot_data$location <- gsub("_", " ", nmds_plot_data$location)
nmds_plot_data$location <- gsub("biopsy", "Biopsy", nmds_plot_data$location)

# Plot the NMDS data
figure_1 <- ggplot(nmds_plot_data, aes(x = NMDS1, y = NMDS2, color = Diagnosis)) +
  geom_point(aes(fill = Diagnosis, shape = Diagnosis),
    size = 1, stroke = 1.2
  ) +
  stat_ellipse(aes(fill = Diagnosis),
    type = "norm",
    level = 0.5,
    alpha = 0.8,
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_text(
    data = permanova_results_wide,
    aes(x = 0, y = Inf, label = label),
    inherit.aes = FALSE,
    color = "black",
    hjust = 0.4,
    vjust = 1.0,
    size = 4,
    fontface = "italic"
  ) +
  labs(x = "NMDS1",
       y = "NMDS2",
       color = "Diagnosis") +
  facet_wrap(~location, ncol = 2, scales = "free") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) + # Add space at top for text
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 8)) + # Add this line
  scale_color_manual(values = c("grey", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  scale_fill_manual(values = c("grey", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  theme_classic() +
  theme(
    axis.line.y = element_line(linetype = 1, linewidth = 1.3, colour = "black"),
    axis.line.x = element_line(linetype = 1, linewidth = 1.3, colour = "black"),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10, color = "black"),
    axis.text.y = element_text(face = "bold", size = 10, color = "black"),
    legend.position = c(0.75, 0.15), # Place legend in the empty bottom-right panel
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(colour = "black", size = 16, face = "bold"),
    legend.text = element_text(colour = "black", size = 14, face = "italic"),
    plot.tag = element_text(face = "bold"),
    strip.background = element_blank(), # Add this line back to remove the box
    strip.text = element_text(face = "bold", size = 12)
  )


print(figure_1)


#define output directory
output_dir <- "/home/maurice/projects/phd/oesphlora/figures/main_figures"

#create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
#define output file
output_file <- file.path(output_dir, "figure_1.pdf")

pdf(output_file, width = 14, height = 8, useDingbats = FALSE)
print(figure_1)
dev.off()