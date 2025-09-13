library(MicrobiomeStat)


asv_table_mat <- as.matrix(asv_table_df)

linda_output <- MicrobiomeStat::linda(
  feature.dat = asv_table_mat,
  meta.dat = metadata_df,
  formula = formula,
  feature.dat.type = "count",
  prev.filter = 0.05,
  mean.abund.filter = 0,
  max.abund.filter = 0,
  is.winsor = TRUE,
  adaptive = TRUE,
  pseudo.cnt = 0.5,
  corr.cut = 0.1,
  n.cores = 1,
  verbose = TRUE
)

linda_results <- as.data.frame(linda_output$output[[1]])
