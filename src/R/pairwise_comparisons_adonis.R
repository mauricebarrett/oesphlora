library(vegan)
library(stats)

# Convert distance matrix to proper format
distance_matrix <- as.dist(distance_matrix)

# Get unique treatment groups
unique_groups <- unique(metadata[[primary_variable]])

print(strata)

# Function to perform pairwise PERMANOVA + Betadisper
pairwise_permanova <- function(group1, group2, metadata, distance_matrix) {
    # Subset metadata for selected groups
    subset_metadata <- metadata[metadata[[primary_variable]] %in% c(group1, group2), ]
    subset_distance <- as.dist(as.matrix(distance_matrix)[rownames(subset_metadata), rownames(subset_metadata)])


    adonis_formula <- as.formula(paste0("subset_distance ~", primary_variable))

    # Perform PERMANOVA
    if (strata == "True") {
    # use patient_id column as strata
    adonis_result <- adonis2(
        adonis_formula,
        data = subset_metadata,
        permutations = 999,
        strata = subset_metadata[["patient_id"]]
    )
    } else {
    # no strata
    adonis_result <- adonis2(
        adonis_formula,
        data = subset_metadata,
        permutations = 999
    )
    }


    # Extract PERMANOVA results
    pseudo_F <- adonis_result$`F`[1]
    p_value <- adonis_result$`Pr(>F)`[1]
    r_squared <- adonis_result$R2[1]

    # Perform betadisper for homogeneity of dispersion
    dispersion_result <- betadisper(subset_distance, subset_metadata[[primary_variable]])
    dispersion_test <- permutest(dispersion_result, permutations = 999)

    # Extract betadisper results
    dispersion_F <- dispersion_test$tab$`F`[1]
    dispersion_p_value <- dispersion_test$tab$`Pr(>F)`[1]

    # Return results
    return(data.frame(
        Group1 = group1, Group2 = group2,
        pseudo_F = pseudo_F, p_value = p_value, r_squared = r_squared,
        dispersion_F = dispersion_F, dispersion_p_value = dispersion_p_value
    ))
}

# Generate all unique pairwise comparisons
pairwise_results <- do.call(rbind, combn(unique_groups, 2, function(pair) {
    pairwise_permanova(pair[1], pair[2], metadata, distance_matrix)
}, simplify = FALSE))

# Convert to data frame
pairwise_results <- as.data.frame(pairwise_results)

# Apply multiple testing correction
pairwise_results$p_value_adjusted <- p.adjust(pairwise_results$p_value, method = "BH")
pairwise_results$dispersion_p_value_adjusted <- p.adjust(pairwise_results$dispersion_p_value, method = "BH")


# Save results to CSV
write.csv(pairwise_results, file = output_file, row.names = TRUE)