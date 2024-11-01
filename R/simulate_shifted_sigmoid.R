simulate_shifted_sigmoid = function(n=4,error_sd=0.05,subject_sd=0,half_temp_grid=FALSE){
  stopifnot(n%%2==0) #we restrict to balanced scenarios only
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

  sample_idx = sample(1:nrow(truepositive_sample_sigmoid_posterior_df),1)
  sample_params = truepositive_sample_sigmoid_posterior_df[sample_idx,] #get posterior params

  abundance_treated = Cetsa(p=sample_params$p_Tmt,k=sample_params$k_Tmt,m=sample_params$m_Tmt,t=temperature)
  abundance_vehicle = Cetsa(p=sample_params$p_Ctrl,k=sample_params$k_Ctrl,m=sample_params$m_Ctrl,t=temperature)

  simdata_wide_treated = MASS::mvrnorm(n=n/2,mu=abundance_treated,Sigma=S)
  colnames(simdata_wide_treated) = NewChannel2Temps_HumanData$Channel


  simdata_wide_vehicle = MASS::mvrnorm(n=n/2,mu=abundance_vehicle,Sigma=S)
  colnames(simdata_wide_vehicle) = NewChannel2Temps_HumanData$Channel

  Subject = paste0("F",seq(1,n,1))

  simdata_wide = rbind(simdata_wide_treated,simdata_wide_vehicle) |> dplyr::as_tibble() |>
   dplyr::mutate(Subject = Subject,Condition=c(rep("treated",n/2),rep("vehicle",n/2)))

  simdata = simdata_wide |> tidyr::pivot_longer(cols=matches("[[:digit:]]+"),names_to="Channel",values_to="Abundance") |>
    dplyr::inner_join(NewChannel2Temps_HumanData,by="Channel") |>
    dplyr::mutate(Group = paste(temperature,Condition,sep="_"),
                  Subject = forcats::as_factor(Subject),
                 Condition = forcats::as_factor(Condition),
                  Condition = forcats::fct_relevel(Condition,"vehicle","treated")) |>
    dplyr::arrange(temperature)

  rho = subject_sd^2/(subject_sd^2 + error_sd^2)
  sample_params = sample_params |> dplyr::select(p_Tmt,k_Tmt,m_Tmt,p_Ctrl,k_Ctrl,m_Ctrl,ATE) |> dplyr::mutate(sigma_e=error_sd,sigma_subject=subject_sd,rho=rho)

  return(list(simdata=simdata,sample_params=sample_params))
}
