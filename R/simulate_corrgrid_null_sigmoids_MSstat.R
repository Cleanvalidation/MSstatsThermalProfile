simulate_corrgrid_null_sigmoids_MSstat = function(runs=1000,n=4,error_sd=0.05,subject_sd=seq(0,0.15,0.025),half_temp_grid=FALSE){
  simdata_list = vector(mode="list",length=length(subject_sd))
  simdata_params_list = vector(mode="list",length=length(subject_sd))

  for(i in 1:length(subject_sd)){
    subSD = subject_sd[i]
    sim = simulate_multiple_null_sigmoids_MSstat(runs=runs,n=n,error_sd=error_sd,subject_sd=subSD,half_temp_grid=half_temp_grid)
    simdata_list[[i]] = sim$simdata
    simdata_params_list[[i]] = sim$simdata_params
    print(paste0("Generated Null Data for Subject SD of ",subSD))
  }

  simdata = dplyr::bind_rows(simdata_list)
  simdata_params = dplyr::bind_rows(simdata_params_list)

  return(list(ProteinLevelData=simdata,params=simdata_params))
}
