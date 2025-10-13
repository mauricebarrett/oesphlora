library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)


primary_variable <- title_list[[1]]
plot_title_names <- title_list[[2]]

daa_df$comparison <- gsub("biopsy_location_", "", daa_df$comparison)


if (primary_variable == "Diagnosis") {
  group_levels <- c("Healthy", "GORD", "BO", "Dysplasia", "OAC", "Metastatic")
} else if (primary_variable == "sample_location") {
  group_levels <- c(
    "1",
    "2",
    "3",
    "4",
    "5"
  )
} else {
  # Error handling for unexpected primary_variable values
  stop(paste("Unexpected primary_variable value:", primary_variable))
}


all_comparisons <- combn(group_levels, 2, function(x) paste(x[1], "vs", x[2]))


all_combos <- expand.grid(
  feature_id = unique(daa_df$feature_id),
  comparison = all_comparisons,
  stringsAsFactors = FALSE
)


# Do a outer join  of plot_df and all_combos to ensure all combinations are present
plot_full_df <- merge(
  daa_df,
  all_combos,
  by = c("comparison", "feature_id"),
  all = TRUE
)

# In log2FoldChange, replace NA values with 0
plot_full_df$log2FoldChange[is.na(plot_full_df$log2FoldChange)] <- 0

# In padj, replace NA values with 1
plot_full_df$padj[is.na(plot_full_df$padj)] <- 1


plot_full_df <- plot_full_df %>%
  mutate(
    star = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  select(feature_id, comparison, log2FoldChange, star, padj)


data_matrix <- dcast(
  plot_full_df,
  feature_id ~ comparison,
  value.var = "log2FoldChange"
)
rownames(data_matrix) <- data_matrix$feature_id
data_matrix <- data_matrix[, -1]

# Replace NA values with 0
data_matrix[is.na(data_matrix)] <- 0

# Cluster features
hc <- hclust(dist(data_matrix, method = "euclidean"), method = "ward.D2")
ordered_features <- rownames(data_matrix)[hc$order]
plot_full_df$feature_id <- factor(
  plot_full_df$feature_id,
  levels = ordered_features
)

# Set comparison as factor to fix x-axis order
plot_full_df$comparison_standard <- factor(
  plot_full_df$comparison,
  levels = all_comparisons
)


plot <- ggplot(plot_full_df, aes(x = comparison_standard, y = feature_id)) +
  geom_tile(aes(fill = log2FoldChange), color = "black") +
  scale_fill_gradient2(
    low = "darkblue",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    space = "Lab",
    na.value = "grey90",
    name = "log2FC"
  ) +
  geom_text(
    aes(label = star),
    color = "black",
    size = 4,
    vjust = 0.7,
    na.rm = TRUE
  ) +
  labs(
    title = plot_title_names,
    x = NULL,
    y = NULL,
    caption = "* <= 0.05, ** <= 0.01, *** <= 0.001"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(colour = "black", size = 10, face = "bold"),
    legend.text = element_text(colour = "black", size = 10, face = "italic"),
    plot.caption = element_text(face = "bold"),
    plot.title = element_text(
      size = 14,
      hjust = 0.5,
      face = "bold",
      color = "black"
    )
  )


# Save the plot
pdf(output_file, width = 14, height = 7.88, useDingbats = FALSE)
print(plot)
dev.off()
