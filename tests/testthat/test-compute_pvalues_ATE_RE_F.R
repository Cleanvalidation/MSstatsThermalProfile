test_that("SCAM splines ATE works", {

  sim=simulate_multiple_shifted_nonsigmoids(runs=10,n=4)
  Sim2MSstat<-Sim2MSstatsTMT(sim)
  scam_ATE=compute_pvalues_ATE_RE_F(Sim2MSstat,workers=2)

  expect_equal(names(scam_ATE),c("Accession","r.sq","dRSS","F_test","F_pvalue","term","contrast","estimate",
               "std.error","statistic","p.value","conf.low","conf.high","F_adjBH","ATE_padjBH"))
})
