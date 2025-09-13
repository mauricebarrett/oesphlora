library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)


# for testing
# data_df <- read.csv("/home/maurice/projects/phd/oesphlora/diversity_metrics/results/diagnosis/Healthy/alpha_diversity_metrics_diagnosis_Healthy.csv")
# create title_list
# title_list <- c("Diagnosis", "Comparison of Alpha Diversity Metrics by Location")

primary_variable <- title_list[[1]]
plot_title_names <- title_list[[2]]

data_df$comparison <- gsub("biopsy_location_", "", data_df$comparison)


data_matrix <- dcast(data_df, metric ~ comparison, value.var = "median_diff")
rownames(data_matrix) <- data_matrix$metric
data_matrix <- data_matrix[, -1]

# Replace NA values with 0
data_matrix[is.na(data_matrix)] <- 0

# Cluster features
hc <- hclust(dist(data_matrix, method = "euclidean"), method = "ward.D2")
ordered_metrics <- rownames(data_matrix)[hc$order]

# Prepare data for plotting
data_df$metric <- factor(data_df$metric, levels = ordered_metrics)


plot <- ggplot(data_df, aes(x = comparison, y = metric)) +
  geom_tile(aes(fill = median_diff), color = "black") + # Moved color outside aes()
  scale_fill_gradient2(
    low = "darkblue", mid = "white", high = "darkred", midpoint = 0,
    space = "Lab", na.value = "grey90", name = "Difference in Medians"
  ) +
  geom_text(aes(label = significance), color = "black", size = 4, vjust = 0.7, na.rm = TRUE, show.legend = FALSE) +
  labs(
    title = plot_title_names,
    x = NULL, y = "Metric",
    caption = "* <= 0.05, ** <= 0.01, *** <= 0.001"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(colour = "black", size = 10, face = "bold"),
    legend.text = element_text(colour = "black", size = 10, face = "italic"),
    plot.caption = element_text(face = "bold"),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold", color = "black")
  )




# Save the plot
pdf(output_file, width = 14, height = 7.88, useDingbats = FALSE)
print(plot)
dev.off()
