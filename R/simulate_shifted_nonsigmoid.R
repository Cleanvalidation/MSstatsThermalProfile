simulate_shifted_nonsigmoid = function(n=4,error_sd=0.05,subject_sd=0,half_temp_grid=FALSE){
  stopifnot(n%%2==0)
  stopifnot(n>=2)

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
  rho = subject_sd^2/(subject_sd^2 + error_sd^2)

  truepos_humandata_proteins = names(nonsigmoid_truepos_humandata_models)
  truepos_zebradata_proteins = names(nonsigmoid_truepos_zebradata_models)

  truepos_proteins = c(truepos_humandata_proteins,truepos_zebradata_proteins)

  #sample protein randomly for now, but can adjust for the distribution of ATEs for overall calculation later

  selected_truepos_protein = sample(truepos_proteins,1,replace=TRUE)

  if(selected_truepos_protein %in% truepos_humandata_proteins){
    selected_truepos_model = nonsigmoid_truepos_humandata_models[[selected_truepos_protein]]
  }else if(selected_truepos_protein %in% truepos_zebradata_proteins){
    selected_truepos_model = nonsigmoid_truepos_zebradata_models[[selected_truepos_protein]]
  }

  treated_temp_df = dplyr::tibble(temperature=temperature,Condition="treated")
  vehicle_temp_df = dplyr::tibble(temperature=temperature,Condition="vehicle")

  abundance_treated = as.vector(predict(selected_truepos_model,newdata=treated_temp_df)) #get mean abundance for treated
  abundance_vehicle = as.vector(predict(selected_truepos_model,newdata=vehicle_temp_df)) #get mean abundance for vehicle


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
                  Condition = forcats::fct_relevel(Condition,"vehicle","treated"),
                  Orig_Protein = selected_truepos_protein,
                  sigma_e=error_sd,sigma_subject=subject_sd,rho=rho) |>
    dplyr::arrange(temperature)

  return(simdata)

}
