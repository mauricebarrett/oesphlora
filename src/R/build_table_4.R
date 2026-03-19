#!/usr/bin/env Rscript

# Build summary table from PERMANOVA results already produced by the pipeline.
# No PERMANOVA is run here; only reads CSVs written by main.py.

# Hardcoded paths (match pipeline output layout)
work_dir <- "/home/mossy/projects/phd/oesphlora"
beta_diagnosis_dir <- file.path(
  work_dir,
  "diversity_metrics",
  "beta_diversity",
  "diagnosis"
)
output_csv <- file.path(work_dir, "table", "beta_diversity_clinical_classification_by_biopsy_location.csv")

read_permanova_values <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  vals <- setNames(as.numeric(df$value), df$metric)
  list(p_value = vals["p_value"], r_squared = vals["r_squared"])
}

metrics <- list(
  list(
    dir_name = "generalized_unifrac",
    file_prefix = "generalized_unifrac",
    col_prefix = "generalized_unifrac"
  ),
  list(
    dir_name = "unweighted_unifrac",
    file_prefix = "unweighted_unifrac",
    col_prefix = "unweighted_unifrac"
  ),
  list(
    dir_name = "jaccard",
    file_prefix = "jaccard",
    col_prefix = "jaccard_index"
  ),
  list(
    dir_name = "rpca",
    file_prefix = "rpca",
    col_prefix = "robust_aitchison"
  ),
  list(
    dir_name = "phylo_rpca",
    file_prefix = "phylo_rpca",
    col_prefix = "phylogenetic_rpca"
  )
)

locations <- paste0("biopsy_location_", 1:5)
rows <- list()

for (loc in locations) {
  loc_num <- sub("biopsy_location_", "", loc)
  row <- list(biopsy_location = loc_num)

  for (m in metrics) {
    path <- file.path(
      beta_diagnosis_dir,
      m$dir_name,
      paste0(m$file_prefix, "_permanova_results_", loc, ".csv")
    )
    stat <- read_permanova_values(path)
    row[[paste0(m$col_prefix, "_p_value")]] <- stat$p_value
    row[[paste0(m$col_prefix, "_r_squared")]] <- round(stat$r_squared, 3)
  }

  rows[[length(rows) + 1]] <- row
}

output_df <- do.call(
  rbind,
  lapply(
    rows,
    function(x) as.data.frame(x, stringsAsFactors = FALSE)
  )
)
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(output_df, output_csv, row.names = FALSE, quote = FALSE)

message("Saved summary table to: ", output_csv)
