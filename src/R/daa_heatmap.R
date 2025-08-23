library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)


primary_variable <- title_list[[1]]
plot_title_names <- title_list[[2]]

group_levels <- c("Healthy", "GORD", "BO", "Dysplasia", "OAC", "Metastatic")


all_comparisons <- combn(group_levels, 2, function(x) paste(x[1], "vs", x[2]))



all_combos <- expand.grid(
  feature_id = unique(daa_df$feature_id),
  comparison = all_comparisons,
  stringsAsFactors = FALSE
)

plot_df <- daa_df %>%
  mutate(
    star = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  select(feature_id, comparison, log2FoldChange, star, padj)


plot_df_full <- left_join(all_combos, plot_df, by = c("feature_id", "comparison"))

# In log2FoldChange, replace NA values with 0
plot_df_full$log2FoldChange[is.na(plot_df_full$log2FoldChange)] <- 0


data_matrix <- dcast(plot_df_full, feature_id ~ comparison, value.var = "log2FoldChange")
rownames(data_matrix) <- data_matrix$feature_id
data_matrix <- data_matrix[, -1]

# Replace NA values with 0
data_matrix[is.na(data_matrix)] <- 0

# Cluster features
hc <- hclust(dist(data_matrix, method = "euclidean"), method = "ward.D2")
ordered_features <- rownames(data_matrix)[hc$order]
plot_df_full$feature_id <- factor(plot_df_full$feature_id, levels = ordered_features)

# Set comparison as factor to fix x-axis order
plot_df_full$comparison_standard <- factor(plot_df_full$comparison, levels = all_comparisons)


plot<-ggplot(plot_df_full, aes(x = comparison_standard, y = feature_id)) +
  geom_tile(aes(fill = log2FoldChange), color = "black") +
  scale_fill_gradient2(
    low = "darkblue", mid = "white", high = "darkred", midpoint = 0,
    space = "Lab", na.value = "grey90", name = "log2FC"
  ) +
  geom_text(aes(label = star), color = "black", size = 4, vjust = 0.7, na.rm = TRUE) +
  labs(
    title = plot_title_names,
    
    x = NULL, y = "Feature",
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

