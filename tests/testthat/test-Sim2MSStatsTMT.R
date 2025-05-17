test_that("Simulations can convert to MSStatsTMT", {
  Sim<-simulate_multiple_null_sigmoids(runs=4)
  Sim<-Sim2MSstatsTMT(Sim$simdata)
  expected_cols <- c("Subject", "Condition", "temperature", "Abundance", "Channel", "Group", "Protein", "BioReplicate")
  missing_cols <- setdiff(expected_cols, names(Sim))

  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }

  Sim <- Sim[, expected_cols]
})
