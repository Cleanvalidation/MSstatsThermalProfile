test_that("Shifted Sigmoid Simulation Correct", {
  set.seed(101)
  sim = simulate_shifted_sigmoid(n=2000,error_sd=1e-8,rho=0)
  true_params = sim$sample_params

  simdata_treated = sim$simdata |> dplyr::filter(Condition=="treated")
  simdata_vehicle = sim$simdata |> dplyr::filter(Condition=="vehicle")

  model_treated = nls(Abundance~Cetsa(p,k,m,temperature),start=list(p=true_params$p_Tmt,k=true_params$k_Tmt,m=true_params$m_Tmt),data=simdata_treated)

  params_treated = coef(model_treated)

  expect_equal(as.numeric(params_treated["p"]),true_params$p_Tmt,tolerance=0.01)
  expect_equal(as.numeric(params_treated["k"]),true_params$k_Tmt,tolerance=0.01)
  expect_equal(as.numeric(params_treated["m"]),true_params$m_Tmt,tolerance=0.01)

  model_vehicle = nls(Abundance~Cetsa(p,k,m,temperature),start=list(p=true_params$p_Ctrl,k=true_params$k_Ctrl,m=true_params$m_Ctrl),data=simdata_vehicle)

  params_vehicle = coef(model_vehicle)

  expect_equal(as.numeric(params_vehicle["p"]),true_params$p_Ctrl,tolerance=0.01)
  expect_equal(as.numeric(params_vehicle["k"]),true_params$k_Ctrl,tolerance=0.01)
  expect_equal(as.numeric(params_vehicle["m"]),true_params$m_Ctrl,tolerance=0.01)


})


test_that("Shifted Sigmoid Simulation Residual Corr Correct", {
  set.seed(101)

  sim = simulate_shifted_sigmoid(n=4000,error_sd=1e-8,rho=0.9)
  true_params = sim$sample_params
  simdata = sim$simdata

  simdata$Tmt = ifelse(simdata$Condition=="treated",1,0)
  simdata$Ctrl = ifelse(simdata$Condition=="vehicle",1,0)

  model = nls(Abundance~Cetsa(p_Tmt,k_Tmt,m_Tmt,temperature)*Tmt +
                Cetsa(p_Ctrl,k_Ctrl,m_Ctrl,temperature)*Ctrl,
              start=list(p_Tmt=true_params$p_Tmt,k_Tmt=true_params$k_Tmt,m_Tmt=true_params$m_Tmt,
                         p_Ctrl=true_params$p_Ctrl,k_Ctrl=true_params$k_Ctrl,m_Ctrl=true_params$m_Ctrl),data=simdata)


  params = coef(model)

  expect_equal(as.numeric(params["p_Tmt"]),true_params$p_Tmt,tolerance=0.01)
  expect_equal(as.numeric(params["k_Tmt"]),true_params$k_Tmt,tolerance=0.01)
  expect_equal(as.numeric(params["m_Tmt"]),true_params$m_Tmt,tolerance=0.01)

  expect_equal(as.numeric(params["p_Ctrl"]),true_params$p_Ctrl,tolerance=0.01)
  expect_equal(as.numeric(params["k_Ctrl"]),true_params$k_Ctrl,tolerance=0.01)
  expect_equal(as.numeric(params["m_Ctrl"]),true_params$m_Ctrl,tolerance=0.01)

  simdata$resids = residuals(model)
  simdata_gdf = simdata |> dplyr::group_split(Subject)
  simdata_corrs = simdata_gdf |>
    purrr::map(.f=function(df) tidyr::expand_grid(res1=df$resids,res2=df$resids) |> dplyr::filter(res1!=res2)) |>
    dplyr::bind_rows()

  rho = cor(simdata_corrs$res1,simdata_corrs$res2)


  expect_equal(rho,0.9,tolerance=0.01)
})
