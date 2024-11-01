test_that("Shifted Non-Sigmoid Simulation Correct", {
  set.seed(99)
  sim = simulate_shifted_nonsigmoid(n=2000,error_sd=1e-15,subject_sd=0)
  selected_protein = unique(sim$Orig_Protein)
  #Check that ATE matches the true for the selected protein
  if(selected_protein %in% names(nonsigmoid_truepos_humandata_models)){
    true_model = nonsigmoid_truepos_humandata_models[[selected_protein]]
  }else if(selected_protein %in% names(nonsigmoid_truepos_zebradata_models)){
    true_model = nonsigmoid_truepos_zebradata_models[[selected_protein]]
  }
  true_marg = marginaleffects::tidy(marginaleffects::comparisons(true_model,var="Condition"))
  true_ATE = true_marg$estimate

  model = scam::scam(Abundance ~ s(temperature,bs="mpd",by=Condition) + Condition,data=sim)
  marg = marginaleffects::tidy(marginaleffects::comparisons(model,var="Condition"))
  ATE_estimate = marg$estimate



  expect_equal(ATE_estimate-true_ATE,0,tolerance=0.01)
})


test_that("Shifted Non-Sigmoid Simulation Residual Corr Correct", {
  set.seed(99)

  #For some reason, cannot set sds too low due to numerical issues making rho to -0.11
  sim = simulate_shifted_nonsigmoid(n=2000,error_sd=0.05,subject_sd=0.15)

  selected_protein = unique(sim$Orig_Protein)
  #Check that ATE matches the true for the selected protein
  if(selected_protein %in% names(nonsigmoid_truepos_humandata_models)){
    true_model = nonsigmoid_truepos_humandata_models[[selected_protein]]
  }else if(selected_protein %in% names(nonsigmoid_truepos_zebradata_models)){
    true_model = nonsigmoid_truepos_zebradata_models[[selected_protein]]
  }
  true_marg = marginaleffects::tidy(marginaleffects::comparisons(true_model,var="Condition"))
  true_ATE = true_marg$estimate

  model = scam::scam(Abundance ~ s(temperature,bs="mpd",by=Condition) + Condition,data=sim)

  marg = marginaleffects::tidy(marginaleffects::comparisons(model,var="Condition"))
  ATE_estimate = marg$estimate

  expect_equal(ATE_estimate-true_ATE,0,tolerance=0.025)

  sim$resids = residuals(model)
  simdata_gdf = sim |> dplyr::group_split(Subject)
  simdata_corrs = simdata_gdf |>
    purrr::map(.f=function(df) tidyr::expand_grid(res1=df$resids,res2=df$resids) |> dplyr::filter(res1!=res2)) |>
    dplyr::bind_rows()

  rho = cor(simdata_corrs$res1,simdata_corrs$res2)

  expect_equal(rho,0.9,tolerance=0.025)
})
