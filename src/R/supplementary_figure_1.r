# Figure 1

library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
suppressPackageStartupMessages(library(ggthemes))

# Load the data
# Load the data
genus_table <- as.matrix(read.csv("/home/maurice/projects/phd/oesphlora/tables/taxa_tables/genus_table.csv", check.names = FALSE, row.names = 1))

# convert to percentage
genus_table_mod <- genus_table[rowSums(genus_table) != 0, ]
genus_table_100 <- apply(genus_table, 2, function(x) { (x / sum(x)) * 100 })




if (nrow(genus_table_100) > 20) {
    len <- nrow(genus_table_100) - 19
    z <- names(sort(rowSums(genus_table_100), decreasing = FALSE)[1:len])
    z <- unique(unlist(z))
    genus_table_filter <- rbind(
        genus_table_100[-which(rownames(genus_table_100) %in% z), , drop = FALSE],
        colSums(genus_table_100[z, , drop = FALSE])
    )
    rownames(genus_table_filter)[nrow(genus_table_filter)] <- "Others"
} else {
    genus_table_filter <- genus_table_100
}

# Convert to long format
genus_table_long <- melt(genus_table_filter , id.vars = "genus", variable.name = "sample_id", value.name = "count")

# Add column names to the genus_table_long dataframe
colnames(genus_table_long) <- c("genus", "sample_id", "count")

# Load metadata
metadata <- read.csv("/home/maurice/projects/phd/oesphlora/metadata/merged_data.csv")

# Merge the data
merged_data <- merge(genus_table_long, metadata, by = "sample_id")

# Aggregate the data by sample_location and genus by finding the mean of the count column
merged_data_aggregated <- merged_data %>%
  group_by(sample_location, genus) %>%
  summarise(count = mean(count), .groups = "drop")


# Add row names as a column manually for both data frames
merged_data_aggregated$sample_location <- factor(merged_data_aggregated$sample_location, levels = c("biopsy_location_1", "biopsy_location_2", "biopsy_location_3", "biopsy_location_4", "biopsy_location_5"))

merged_data_aggregated$sample_location <- gsub("_", " ", merged_data_aggregated$sample_location)
merged_data_aggregated$sample_location <- gsub("biopsy", "Biopsy", merged_data_aggregated$sample_location)

# Reorder variable (asv) column based on value (Abundance in percentage)
merged_data_aggregated$genus <- reorder(merged_data_aggregated$genus, merged_data_aggregated$count)

# Plot genus level data
taxa_bar_plot <- ggplot(data = merged_data_aggregated, aes(x = sample_location, y = count, fill = genus)) +
            geom_bar(stat = "identity", position = "stack", lwd = 0) +
            scale_y_continuous(breaks = seq(0, 100, 10), expand = expansion(mult = c(0, 0))) +
            scale_fill_manual(values = tableau_color_pal(palette = "Tableau 20", direction = -1)(20)) +
            labs(title = NULL,
                 x = "", y = "Relative abundance (%)") +
            guides(fill = guide_legend(title = "Genus", reverse = FALSE)) +
            theme(panel.background = element_rect(fill = 'white'),
                  plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(face = "bold", size = 12, hjust = 0.5, family = "Helvetica"),
                  axis.line.y = element_line(linetype = 1, linewidth = 0.5, colour = 'black'),
                  axis.line.x = element_line(linetype = 1, linewidth = 0.5, colour = 'black'),
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "bold", size = 12),
                  axis.text.y = element_text(face = "bold", size = 12, color = "black"),
                  axis.title.y = element_text(face = "bold", size = 12),
                  axis.title.x = element_text(face = "bold", size = 12),
                  legend.position = "right",
                  legend.title = element_text(colour = "black", size = 12, face = "bold"),
                  legend.text = element_text(colour = "black", size = 10, face = "italic"))


print(taxa_bar_plot)


#define output directory
output_dir <- "/home/maurice/projects/phd/oesphlora/figures/supplementary_figures"

#create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
#define output file
output_file <- file.path(output_dir, "supplementary_figure_1.pdf")

pdf(output_file, width = 14, height = 8, useDingbats = FALSE)
print(taxa_bar_plot)
dev.off()