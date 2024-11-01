test_that("Simulations can convert to MSStatsTMT", {
  Sim<-simulate_multiple_null_sigmoids(runs=4)
  Sim<-Sim2MSstatsTMT(Sim$simdata)
  testthat::expect_equal(names(Sim),c("Subject",
                            "Condition",
                            "temperature",
                            "Abundance",
                            "Channel",
                            "Group",
                            "Protein",
                            "BioReplicate",
                            "Mixture",
                            "Run",
                            "Experiment",
                            "Fraction",
                            "TechRepMixture"))
})
