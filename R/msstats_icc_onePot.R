msstats_icc_onePot<-function(summarised.proteins,contrast.matrix){
  if (!requireNamespace("MSstatsTMT", quietly = TRUE)) {
    stop("The MSstatsTMT package is required but not installed.")
  }
  if (!requireNamespace("MSstats", quietly = TRUE)) {
    stop("The MSstats package is required but not installed.")
  }
  if (!requireNamespace("MSstatsConvert", quietly = TRUE)) {
    stop("The MSstatsConvert package is required but not installed.")
  }

  #define the contrast matrix
  comparison<-contrast.matrix
  #use MSstatsTMT statistical models to evaluate contrasts at the specific temperatures
  DIM_MSstats<-MSstatsTMT::groupComparisonTMT(
    summarised.proteins,
    contrast.matrix = comparison,
    moderated = FALSE,
    adj.method = "BH",
    remove_norm_channel = TRUE,
    remove_empty_channel = FALSE,
    save_fitted_models = TRUE,
    use_log_file = FALSE,
    append = FALSE,
    verbose = TRUE,
    log_file_path = NULL
  )
  #this calculates the residuals from MSstats output if the class is lmm
  VarCor_MSstats<-lapply(DIM_MSstats$FittedModel,function(x) tryCatch(print(lme4::VarCorr(x,comp="Variance"))|>as.data.frame(),error=function(x){return(NA)}))
  #Check if all models are LMM and show how many proteins converged
  VarCor_fit<-VarCor_MSstats|>purrr::keep(function(x) !is.null(ncol(x)))
  print(paste0(length(VarCor_fit)," out of ",length(VarCor_MSstats)," proteins converged with MSstatsTMT"))
  #transform the original data to a list to append the variance components
  protein<-data.frame(Protein=names(VarCor_fit))
  proteinList<-summarised.proteins$ProteinLevelData[summarised.proteins$ProteinLevelData$Protein %in% protein$Protein,]
  proteinList<-proteinList|>dplyr::group_by(Protein)|>dplyr::group_split()

  #append to the protein list the variance components and icc values
  out_df<-purrr::map2(proteinList,VarCor_fit,
                      function(x,y)x|>dplyr::mutate(
                        sigma_bio=y$sdcor[1],
                        sigma_re=y$sdcor[2],
                        icc=y$vcov[1]/(y$vcov[1]+y$vcov[2]
                        )))

  icc_data<-lapply(out_df,function(x) x[1,])|>dplyr::bind_rows()#save the first row for every protein to do a histogram of icc values
  #QC histogram for ICC values
  hist(icc_data$icc,20,main="TPP proc. Dataset   (all proteins)",xlab="ICC",col="#D95F0E")
  #add a treatment column
  out_df<-dplyr::bind_rows(out_df)|>dplyr::inner_join(summarised.proteins$ProteinLevelData|>dplyr::select(Run,Condition)|>dplyr::distinct())
  hist(icc_data$sigma_bio,20,main="TPP proc. Dataset  (all proteins)",xlab="sigma_bio",col="#D95F0E")

  hist(icc_data$sigma_re,20,main="TPP proc. Dataset  (all proteins)",xlab="sigma_tech",col="#D95F0E")
  summary(icc_data$icc)
  summary(icc_data$sigma_bio)

  summary(icc_data$sigma_tech)
  df_with_variance_fittedModel <- list("df_with_variance" = out_df, "fitted_model" = DIM_MSstats)
  return (df_with_variance_fittedModel)
}
