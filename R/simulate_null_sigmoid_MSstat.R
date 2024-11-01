
simulate_null_sigmoid_MSstat = function(n=4,error_sd=0.05,subject_sd=0,half_temp_grid=FALSE){
  stopifnot(n%%2==0)
  temperature = c(37,41,44,47,50,53,56,59,63,67)
  NewChannel2Temps_HumanData = Channel2Temps_HumanData[Channel2Temps_HumanData$temperature %in% temperature,]
  if(half_temp_grid==TRUE){
    temperature = rep(temperature[c(1,3,5,7,9)],2)
    NewChannel2Temps_HumanData$temperature = temperature
  }
  Ntemps = length(temperature)

  # random effects Z is a column vec of 1s for RE intercept, ZZ^T is an Ntemps x Ntemps matrix of 1s where k = # of temps
  # G is just RE intercept variance (aka subject_sd^2)--scalar because only RE intercept in this case

  ZGZT = subject_sd^2 * matrix(rep(1,Ntemps^2),nrow=Ntemps,ncol=Ntemps)

  S = ZGZT + error_sd^2 * diag(Ntemps)

  sample_idx = sample(1:nrow(null_sigmoid_posterior_df_MSstats),1)
  sample_params = null_sigmoid_posterior_df_MSstats[sample_idx,] #get posterior params

  abundance = Cetsa(sample_params$p,sample_params$k,sample_params$m,t=temperature)

  simdata_wide = MASS::mvrnorm(n=n,mu=abundance,Sigma=S)
  colnames(simdata_wide) = NewChannel2Temps_HumanData$Channel
  Subject = paste0("F",seq(1,n,1))
  simdata_wide = dplyr::as_tibble(simdata_wide) |> dplyr::mutate(Subject = Subject,Condition=c(rep("treated",n/2),rep("vehicle",n/2)))

  simdata = simdata_wide |> tidyr::pivot_longer(cols=matches("[[:digit:]]+"),names_to="Channel",values_to="Abundance") |>
    dplyr::inner_join(NewChannel2Temps_HumanData,by="Channel") |>
    dplyr::mutate(Group = paste(temperature,Condition,sep="_"),
                  Subject = forcats::as_factor(Subject),
                  Condition = forcats::as_factor(Condition),
                  Condition = forcats::fct_relevel(Condition,"vehicle","treated")) |>
    dplyr::arrange(temperature)

  rho = subject_sd^2/(subject_sd^2 + error_sd^2)
  sample_params = sample_params |> dplyr::select(p,k,m) |> dplyr::mutate(sigma_e=error_sd,sigma_subject=subject_sd,rho=rho)


  return(list(simdata=simdata,sample_params=sample_params))


  return(simdata)
}
