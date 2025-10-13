library(vegan)

set.seed(42) # Set seed for reproducibility

distance_matrix <- as.dist(distance_matrix)

# Dynamically create the formula using primary_variable
adonis_formula <- as.formula(paste("distance_matrix ~", primary_variable))


if (strata == "True") {
  # use patient_id column as strata
  adonis_result <- adonis2(
    adonis_formula,
    data = metadata,
    permutations = 999,
    strata = metadata[["patient_id"]]
  )
} else {
  # no strata
  adonis_result <- adonis2(
    adonis_formula,
    data = metadata,
    permutations = 999
  )
}

# Extract specific PERMANOVA results
pseudo_F <- adonis_result$F[1]
p_value <- adonis_result$`Pr(>F)`[1]
r_squared <- adonis_result$R2[1]

# Dynamically select the grouping column from metadata
grouping_vector <- metadata[[primary_variable]]

# Perform betadisper
dispersion_result <- betadisper(distance_matrix, grouping_vector)
dispersion_test <- permutest(dispersion_result, permutations = 999)

# Extract betadisper results
dispersion_F <- dispersion_test$tab$F[1]
dispersion_p_value <- dispersion_test$tab$`Pr(>F)`[1]

results_df <- data.frame(
  metric = c(
    "pseudo_F",
    "p_value",
    "r_squared",
    "dispersion_F",
    "dispersion_p_value"
  ),
  value = c(pseudo_F, p_value, r_squared, dispersion_F, dispersion_p_value)
)
