test_that("Sigmoid Simulation Correct", {
  set.seed(101)
  sim = simulate_null_sigmoid(n=1000,error_sd=1e-8,rho=0)
  true_params = sim$sample_params
  simdata = sim$simdata

  model = nls(Abundance~Cetsa(p,k,m,temperature),start=list(p=true_params$p,k=true_params$k,m=true_params$m),data=simdata)

  params = coef(model)
  expect_equal(as.numeric(params["p"]),true_params$p,tolerance=0.01)
  expect_equal(as.numeric(params["k"]),true_params$k,tolerance=0.01)
  expect_equal(as.numeric(params["m"]),true_params$m,tolerance=0.01)
})

test_that("Sigmoid Simulation Residual Corr Correct", {
  set.seed(101)

  sim = simulate_null_sigmoid(n=1000,error_sd=1e-8,rho=0.9)
  true_params = sim$sample_params
  simdata = sim$simdata

  model = nls(Abundance~Cetsa(p,k,m,temperature),start=list(p=true_params$p,k=true_params$k,m=true_params$m),data=simdata)

  params = coef(model)

  expect_equal(as.numeric(params["p"]),true_params$p,tolerance=0.01)
  expect_equal(as.numeric(params["k"]),true_params$k,tolerance=0.01)
  expect_equal(as.numeric(params["m"]),true_params$m,tolerance=0.01)

  simdata$resids = residuals(model)
  simdata_gdf = simdata |> dplyr::group_split(Subject)
  simdata_corrs = simdata_gdf |>
    purrr::map(.f=function(df) tidyr::expand_grid(res1=df$resids,res2=df$resids) |> dplyr::filter(res1!=res2)) |>
    dplyr::bind_rows()

  rho = cor(simdata_corrs$res1,simdata_corrs$res2)

  expect_equal(rho,0.9,tolerance=0.01)
})
