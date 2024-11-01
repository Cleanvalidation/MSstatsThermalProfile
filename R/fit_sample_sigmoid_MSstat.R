fit_sample_sigmoid_MSstat = function(){
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

  sigmoid_stanmod = rstan::stan_model("R/sigmoid_treatment_model.stan")
  print("Loaded Sample Sigmoid Model, beginning fitting to sample Q6NV46 Zebrafish data MSstatsTMT processed")


  Q6NV46_sampledata_MSstat = list(N=nrow(Zebra_Q6NV46_MSstat),Abundance=Zebra_Q6NV46_MSstat$Abundance,
                                  Temperature=Zebra_Q6NV46_MSstat$temperature,
                                  Treatment=ifelse(Zebra_Q6NV46_MSstat$Treatment=="treated",1,0))

  results = rstan::sampling(sigmoid_stanmod,data=Q6NV46_sampledata_MSstat,control=list(adapt_delta=0.999,max_treedepth=15,stepsize=0.01),iter=2500,
                            chains=4,cores=4,seed=101)

  posterior_df = posterior::as_draws_df(results) |> dplyr::as_tibble()

  return(list(model=results,posterior_df = posterior_df))

}
