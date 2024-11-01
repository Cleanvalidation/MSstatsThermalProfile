fit_sample_null_sigmoid = function(){
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop(
      "Package \"rstan\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop(
      "Package \"posterior\" must be installed to use this function.",
      call. = FALSE
    )
  }

  sigmoid_stanmod = rstan::stan_model("R/sigmoid_null_model.stan")
  print("Loaded Sample Sigmoid Model, beginning fitting to sample O00267 Human data")


  O00267_sampledata = list(N=nrow(Human_O00267),Abundance = Human_O00267$Abundance,Temperature=Human_O00267$temperature,
                           Treatment=ifelse(stringr::str_detect(Human_O00267$Condition,"treated"),1,0),
                           NSubjects = length(unique(Human_O00267$Mixture)),
                           Subject_ID = as.numeric(stringr::str_remove(Human_O00267$Mixture,"F")))

  results = rstan::sampling(sigmoid_stanmod,data=O00267_sampledata,control=list(adapt_delta=0.999,max_treedepth=15,stepsize=0.01),iter=2500,
                            chains=4,cores=4,seed=101)

  posterior_df = posterior::as_draws_df(results) |> dplyr::as_tibble()

  return(list(model=results,posterior_df = posterior_df))

}
