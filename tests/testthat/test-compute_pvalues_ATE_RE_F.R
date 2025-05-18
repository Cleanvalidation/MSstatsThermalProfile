test_that("SCAM splines ATE works", {

  sim=simulate_multiple_shifted_nonsigmoids(runs=10,n=4)
  Sim2MSstat<-Sim2MSstatsTMT(sim) |>
    dplyr::mutate(Accession=Protein)
  scam_ATE=compute_pvalues_ATE_RE_F(Sim2MSstat,workers=2)

  required_cols <- c("r.sq", "dRSS", "F_test", "F_pvalue", "term", "contrast",
                     "estimate", "std.error", "statistic", "p.value", "conf.low",
                     "conf.high", "F_adjBH", "ATE_padjBH")
  expect_true(all(required_cols %in% names(scam_ATE)))

})
