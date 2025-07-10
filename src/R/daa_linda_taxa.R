library(MicrobiomeStat)

linda_output = MicrobiomeStat::linda(
    taxon_table_df,
    metadata_df,
    formula = formula,
    feature.dat.type = 'count',
    prev.filter = 0.10,
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