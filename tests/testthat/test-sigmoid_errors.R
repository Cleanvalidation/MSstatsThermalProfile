test_that("sigmoid error partitioning function works", {
  data("TPPnorm_Human_Proteins",package="MSstatsThermalProfiler")
  all=TPPnorm_Human_Proteins

  errors <-sigmoid_errors(all,proteinName=c("P36507"),cond="vehicle",n_protein=1000)
})
