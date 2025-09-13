################################################################################
# MICROBIAL CO-OCCURRENCE NETWORK ANALYSIS USING SPARCC
################################################################################
#
# Description:
#   This R script performs microbial co-occurrence network analysis using SparCC
#   (Sparse Correlations for Compositional data) correlation estimation. It
#   creates correlation matrices, hierarchical clustering heatmaps, and network
#   visualizations with community detection and hub identification.
#
# Input Requirements (passed from Python via rpy2):
#   - taxa_table_matrix: Filtered taxa abundance matrix (taxa x samples)
#   - output_dir: Output directory path for saving plots and results
#
# Optional Parameters:
#   - correlation_threshold: Minimum correlation threshold (default: 0.1)
#   - nboot: Number of bootstrap iterations (default: 100)
#   - alpha: Significance level for edge selection (default: 0.05)
#   - centrality_percentile: Percentile threshold for hub identification (default: 0.95)
#
# Output Files:
#   - sparcc_weighted_igraph_network_plot_with_netcomi.pdf: Network visualization
#   - correlation_matrix.csv: SparCC correlation matrix (if saved)
#   - network_stats.txt: Network summary statistics (if saved)
#
# Dependencies:
#   - ggplot2: Plotting framework
#   - SpiecEasi: Alternative correlation methods (if needed)
#   - Matrix: Sparse matrix operations
#   - igraph: Graph objects and network metrics
#   - ggraph: Network visualization with ggplot2
#   - reshape2: Data manipulation
#   - Cairo: PDF output with better quality
#
# Analysis Steps:
#   1. Data preprocessing and matrix transformation
#   2. SparCC correlation estimation with bootstrap validation
#   3. Hierarchical clustering of correlation matrix
#   4. Network construction with threshold filtering
#   5. Community detection using fast greedy algorithm
#   6. Eigenvector centrality calculation for hub identification
#   7. Network visualization with nodes colored by community membership
#   8. Export high-quality PDF plot
#
# Note: This script is designed to be called from Python using rpy2 and expects
#       specific variable names to be available in the R environment.
#
################################################################################



# Load libraries
library(ggplot2)
library(SpiecEasi)
library(igraph)
library(ggraph)
library(reshape2)
library(Cairo)
library(ggdendro)
library(cowplot)

primary_variable <- title_list[[1]]
location <- title_list[[2]]
rank <- title_list[[3]]


# Convert dataframe to matrix
taxa_table_matrix <- as.matrix(taxa_table_df)


# Define path to out file for correlation matrix
correlation_matrix_file <- file.path(output_dir, paste0(rank, "_sparcc_correlation_matrix_location_", location, ".csv"))

# Define path to file for p-value matrix
p_value_matrix_file <- file.path(output_dir, paste0(rank, "_sparcc_p_value_matrix_location_", location, ".csv"))

# Define path to output file for heatmap plot
heatmap_plot_file <- file.path(output_dir, paste0(rank, "_sparcc_heatmap_location_", location, ".pdf"))

# Define path to output file for network plot
network_plot_file <- file.path(output_dir, paste0(rank, "_sparcc_network_plot_location_", location, ".pdf"))


# transpose the matrix
taxa_table_matrix <- t(taxa_table_matrix)

# Ensure the matrix is numeric
mode(taxa_table_matrix) <- "numeric"

# print dimensions of the matrix
print(paste0("Dimensions of taxa_table_matrix: ", dim(taxa_table_matrix)[1], " rows, ", dim(taxa_table_matrix)[2], " columns."))

print(paste0("Running SparCC on ", rank, " data for location: ", location, "..."))

sparcc_result <- sparcc(taxa_table_matrix, iter = 20, inner_iter = 10, th = 0.1)


print(paste0("SparCC analysis completed for ", rank, " data for location: ", location, "."))

correlation_matrix <- as.matrix(sparcc_result$Cor)


feature_names <- colnames(taxa_table_matrix)
rownames(correlation_matrix) <- feature_names
colnames(correlation_matrix) <- feature_names


# Faster parameters for development/testing
boot <- 50
ncpus <- 16
sparcc_params <- list(
  iter = 20, # Reduced from 50 iterations
  inner_iter = 10, # Reduced from 20 inner iterations
  th = 0.1 # Keep threshold the same
)

print(paste0("Running SparCC bootstrap with ", boot, " iterations using ", ncpus, " cores..."))

# Perform SparCC bootstrap analysis
sb <- sparccboot(
  data = taxa_table_matrix,
  sparcc.params = sparcc_params,
  R = boot,
  ncpus = ncpus,
)


print(paste0("SparCC bootstrap completed for ", rank, " data for location: ", location, "."))


# Calculate p-values for the SparCC bootstrap results
sparccboot_result <- pval.sparccboot(sb, sided = "both")

# Extract p-values and adjust them using Benjamini-Hochberg method
p_adj_vec <- p.adjust(sparccboot_result$pvals, method = "BH")

# Create a matrix from the adjusted p-values
p_value_matrix <- triu2diag(p_adj_vec, diagval = 0)


# Set row and column names for the p-value matrix
rownames(p_value_matrix) <- feature_names
colnames(p_value_matrix) <- feature_names


# Save the correlation matrix to a CSV file
write.csv(correlation_matrix, file = correlation_matrix_file, row.names = TRUE)


# Save the p-value matrix to a CSV file
write.csv(p_value_matrix, file = p_value_matrix_file, row.names = TRUE)


# # For the second graph, use the SAME seeds:
# net_sparcc <- netConstruct(
#   data       = mat,
#   dataType   = "counts",
#   measure    = "sparcc",
#     measurePar  = list(
#     iter       = 50,
#     inner_iter = 20,
#     th         = 0.1
#   ),
#   # SparCC has its own compositional/zero handling — don’t double-normalize
#   normMethod = "none",
#   zeroMethod = "none",
#   sparsMethod = "bootstrap",    # Friedman & Alm-style edge selection
#   nboot       = 50,            # increase (e.g., 500–1000) for publication runs
#   adjust      = "adaptBH",      # multiple-testing control for boot p-values
#   alpha       = 0.05,
#   # light filters; tweak to your data size
#   # filtSamp   = "totalReads",  filtSampPar = list(totalReads = 1000),
#   # filtTax    = "numbSamp",    filtTaxPar  = list(numbSamp = 5),
#   verbose = 2, seed=123
# )


###################################
# Heatmap of the correlation matrix
###################################




heatmap_data <- correlation_matrix

# Remove pairs with p-value > 0.05
heatmap_data[p_value_matrix > 0.05] <- 0

# Perform hierarchical clustering on rows and columns
row_dist <- dist(correlation_matrix, method = "euclidean")
col_dist <- dist(t(correlation_matrix), method = "euclidean")

row_clust <- hclust(row_dist, method = "complete")
col_clust <- hclust(col_dist, method = "complete")

# Reorder the matrix based on clustering FIRST
correlation_matrix_clustered <- correlation_matrix[row_clust$order, col_clust$order]

# Convert clustered matrix to long format for ggplot BEFORE using it
correlation_melted <- melt(correlation_matrix_clustered)
names(correlation_melted) <- c("feature1", "feature2", "Correlation")

# Convert factor levels to maintain clustering order
correlation_melted$feature1 <- factor(correlation_melted$feature1,
  levels = rownames(correlation_matrix_clustered)
)
correlation_melted$feature2 <- factor(correlation_melted$feature2,
  levels = colnames(correlation_matrix_clustered)
)

# Create dendrograms as grobs
row_dendro <- as.dendrogram(row_clust)
col_dendro <- as.dendrogram(col_clust)

# Convert dendrograms to ggplot objects
row_dendro_plot <- ggdendrogram(rev(row_dendro), rotate = TRUE, size = 0.5) +
  scale_y_reverse() + # Flip the visual direction
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


# NOW create the main heatmap using the properly defined correlation_melted
heatmap_main <- ggplot(correlation_melted, aes(x = feature1, y = feature2, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "red", mid = "white", high = "blue",
    midpoint = 0,
    name = "SparCC\nCorrelation",
    limits = c(-1, 1)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    legend.title = element_text(size = 10, face = "bold"),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )




# Alternative: More precise control with cowplot - ROW DENDROGRAM ON LEFT
final_plot <- ggdraw() +
  draw_plot(row_dendro_plot, x = 0, y = 0.07, width = 0.2, height = 0.75) + # Negative y value
  draw_plot(heatmap_main, x = 0.2, y = 0, width = 0.8, height = 0.9) +
  draw_label(paste0("SparCC Correlation Heatmap - ", rank, " (", location, ")"),
    x = 0.6, y = 0.95, fontface = "bold", size = 14
  )


# Save the precisely aligned version
save_plot(heatmap_plot_file,
  final_plot,
  base_width = 16,
  base_height = 12
)

plot_grid(row_dendro_plot, heatmap_main, align = "h", rel_widths = c(1, 1))


##################
# Network analysis
##################

# Set any values below the threshold to zero
correlation_matrix[abs(correlation_matrix) < 0.2] <- 0

# # Filter the correlation matrix based on the adjusted p-values
correlation_matrix[p_value_matrix > 0.05] <- 0

diag(correlation_matrix) <- 0

# Ensure the matrix is numeric
mode(correlation_matrix) <- "numeric"

# Create igraph object directly from matrix with names preserved
g <- graph_from_adjacency_matrix(correlation_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)

# Extract the largest connected component (LCC)
components <- components(g)
largest_component_id <- which.max(components$csize)
vertices_in_lcc <- which(components$membership == largest_component_id)

# Create subgraph with only the largest connected component
g_lcc <- induced_subgraph(g, vids = vertices_in_lcc)

# Continue analysis with the LCC only
g <- g_lcc # Replace g with the LCC for all subsequent analyses

# Create new edge attributes called corr
E(g)$corr <- E(g)$weight
# Create new edge attributes called abs_corr
E(g)$abs_corr <- abs(E(g)$corr)
# Replace NA values with 0
E(g)$corr[is.na(E(g)$corr)] <- 0
# Replace NA values with 0 in abs_corr
E(g)$abs_corr[is.na(E(g)$abs_corr)] <- 0

# Use cluster_fast_greedy to detect communities
com <- cluster_fast_greedy(g, weights = E(g)$abs_corr)

# Assign community membership as a vertex attribute
V(g)$community <- factor(membership(com))

# Calculate eigenvector centrality for the LCC
vertex_evc_scores <- eigen_centrality(g, directed = FALSE, weights = E(g)$abs_corr)$vector
vertex_evc_scores[is.na(vertex_evc_scores)] <- 0

# Optional normalization
if (max(vertex_evc_scores) > 0) {
  vertex_evc_scores <- vertex_evc_scores / max(vertex_evc_scores)
}

# Assign as a vertex attribute for plotting/analysis
V(g)$eig <- vertex_evc_scores

# Identify hub nodes based on 95th percentile of eigenvector centrality
centrality_threshold <- quantile(vertex_evc_scores, 0.95) # Top 5% as hubs
V(g)$is_hub <- vertex_evc_scores > centrality_threshold

# Keep the signed values in E(g)$corr for coloring, but use abs weights for layout algos
E(g)$weight <- E(g)$abs_corr

# --- Plot with ggraph ---
set.seed(123)
p1 <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = abs_corr, color = corr > 0), alpha = 0.8) +
  scale_edge_width(range = c(0.5, 3), name = "Correlation |r|") +
  scale_edge_color_manual(
    values = c(`TRUE` = "#009900", `FALSE` = "red"),
    labels = c(`TRUE` = "Positive", `FALSE` = "Negative"),
    name   = "Estimated association"
  ) +
  geom_node_point(aes(size = eig, color = community, shape = is_hub), alpha = 0.9) +
  scale_size_continuous(range = c(3, 12), guide = "none") +
  scale_color_brewer(type = "qual", palette = "Set3", name = "Community", guide = "none") +
  scale_shape_manual(
    values = c(`TRUE` = 15, `FALSE` = 19),
    name = "Node Type",
    labels = c(`TRUE` = "Hub", `FALSE` = "Regular")
  ) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3, fontface = "bold") +
  theme_graph() +
  labs(
    title = "Co-occurrence Network analysis of SparCC associations (Largest Connected Component)",
    subtitle = paste0(rank, " - ", vcount(g), " nodes, ", ecount(g), " edges")
  ) +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 12, hjust = 0.5, face = "italic"),
    legend.position = "left",
    legend.title = element_text(colour = "black", size = 10, face = "bold"),
    legend.text = element_text(colour = "black", size = 8)
  )

# --- Save ---
CairoPDF(file = network_plot_file, width = 14, height = 8)
print(p1)
dev.off()
