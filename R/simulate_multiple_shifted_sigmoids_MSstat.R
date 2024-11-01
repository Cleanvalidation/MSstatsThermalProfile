simulate_multiple_shifted_sigmoids_MSstat = function(runs=1000,n=4,error_sd=0.05,subject_sd=0,half_temp_grid=FALSE){
  simdata_list = vector(mode="list",length=runs)
  simdata_params_list = vector(mode="list",length=runs)

  for(i in 1:runs){
    Protein = paste("Shifted_SubSD",subject_sd,"Sim",i,sep="_")

    sim = simulate_shifted_sigmoid_MSstat(n=n,error_sd=error_sd,subject_sd = subject_sd,half_temp_grid=half_temp_grid)
    simdata_list[[i]] = sim$simdata
    simdata_list[[i]]$Protein = Protein

    simdata_params_list[[i]] = sim$sample_params
    simdata_params_list[[i]]$Protein = Protein
  }

  simdata = dplyr::bind_rows(simdata_list)
  simdata_params = dplyr::bind_rows(simdata_params_list)

  return(list(simdata=simdata,simdata_params=simdata_params))
}
