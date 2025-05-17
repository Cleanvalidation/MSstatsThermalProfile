#' @importFrom scam predict


simulate_null_nonsigmoid_MSstat = function(n=4,error_sd=0.05,subject_sd=0,half_temp_grid=FALSE){
  stopifnot(n%%2==0)
  stopifnot(n>=2)
  temperature = c(37,41,44,47,50,53,56,59,63,67)
  NewChannel2Temps_HumanData = Channel2Temps_HumanData[Channel2Temps_HumanData$temperature %in% temperature,]
  if(half_temp_grid==TRUE){
    temperature = rep(temperature[c(1,3,5,7,9)],2)
    NewChannel2Temps_HumanData$temperature = temperature
  }
  Ntemps = length(temperature)

  ZGZT = subject_sd^2 * matrix(rep(1,Ntemps^2),nrow=Ntemps,ncol=Ntemps)

  S = ZGZT + error_sd^2 * diag(Ntemps)
  rho = subject_sd^2/(subject_sd^2 + error_sd^2)

  true_neg_proteins = names(nonsigmoid_trueneg_humandata_models_MSstats)

  selected_true_neg_protein = sample(true_neg_proteins,1,replace=TRUE) #sample a protein at random

  selected_true_neg_model = nonsigmoid_trueneg_models[[selected_true_neg_protein]] #get model corresponding to that protein

  temp_df = dplyr::tibble(temperature=temperature)

  abundance = as.vector(predict(selected_true_neg_model,newdata=temp_df)) #get mean abundance

  simdata_wide = MASS::mvrnorm(n=n,mu=abundance,Sigma=S)
  colnames(simdata_wide) = NewChannel2Temps_HumanData$Channel

  Subject = paste0("F",seq(1,n,1))
  simdata_wide = dplyr::as_tibble(simdata_wide) |> dplyr::mutate(Subject = Subject,Condition=c(rep("treated",n/2),rep("vehicle",n/2)))

  simdata = simdata_wide |> tidyr::pivot_longer(cols=matches("[[:digit:]]+"),names_to="Channel",values_to="Abundance") |>
    dplyr::inner_join(NewChannel2Temps_HumanData,by="Channel") |>
    dplyr::mutate(Group = paste(temperature,Condition,sep="_"),
                  Subject = forcats::as_factor(Subject),
                  Condition = forcats::as_factor(Condition),
                  Condition = forcats::fct_relevel(Condition,"vehicle","treated"),
                  Orig_Protein = selected_true_neg_protein,
                  sigma_e=error_sd,sigma_subject=subject_sd,rho=rho) |>
    dplyr::arrange(temperature)

  return(simdata)

}
