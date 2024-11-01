test_that("Non-Sigmoid Simulation Correct", {
  set.seed(101)
  sim = simulate_null_nonsigmoid(n=2000,error_sd=1e-15,subject_sd=0)

  #Check that ATE is 0
  model = scam::scam(Abundance ~ s(temperature,bs="mpd",by=Condition) + Condition,data=sim)
  marg = marginaleffects::tidy(marginaleffects::comparisons(model,var="Condition"))
  ATE_estimate = marg$estimate

  expect_equal(ATE_estimate,0)
})


test_that("Non-Sigmoid Simulation Residual Corr Correct", {
  set.seed(101)

  #For some reason, cannot set sds too low due to numerical issues making rho to -0.11
  sim = simulate_null_nonsigmoid(n=1000,error_sd=0.05,subject_sd=0.15)

  model = scam::scam(Abundance ~ s(temperature,bs="mpd",by=Condition) + Condition,data=sim)

  marg = marginaleffects::tidy(marginaleffects::comparisons(model,var="Condition"))
  ATE_estimate = marg$estimate

  expect_equal(ATE_estimate,0,tolerance=0.03)

  sim$resids = residuals(model)
  simdata_gdf = sim |> dplyr::group_split(Subject)
  simdata_corrs = simdata_gdf |>
    purrr::map(.f=function(df) tidyr::expand_grid(res1=df$resids,res2=df$resids) |> dplyr::filter(res1!=res2)) |>
    dplyr::bind_rows()

  rho = cor(simdata_corrs$res1,simdata_corrs$res2)

  expect_equal(rho,0.9,tolerance=0.01)
})

