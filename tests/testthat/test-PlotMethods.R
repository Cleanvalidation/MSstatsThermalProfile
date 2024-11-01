test_that("dataset has the correct format for plots", {
  sim = simulate_null_sigmoid(n=1000,error_sd=1e-8,subject_sd=0)
  df<-sim$simdata
  ch_names<-names(MSstatsThermalProfiler::Zebra_Q6NV46|>dplyr::rename(Subject="BioReplicate")|>
                    dplyr::mutate(Group=stringr::str_extract(Condition,"[[:lower:]]+"))|>
                    dplyr::select(Subject,Condition,Channel,Abundance,temperature,Group))

  expect_equal(names(df),ch_names)
})
