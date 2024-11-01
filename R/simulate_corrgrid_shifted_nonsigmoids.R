simulate_corrgrid_shifted_nonsigmoids = function(runs=1000,n=4,error_sd=0.05,subject_sd=seq(0,0.15,0.025),half_temp_grid=FALSE){
  simdata_list = vector(mode="list",length=length(subject_sd))

  for(i in 1:length(subject_sd)){
    subSD = subject_sd[i]
    sim = simulate_multiple_shifted_nonsigmoids(runs=runs,n=n,error_sd=error_sd,subject_sd=subSD,half_temp_grid=half_temp_grid)
    simdata_list[[i]] = sim
    print(paste0("Generated Shifted Data for Subject SD of ",subSD))
  }

  simdata = dplyr::bind_rows(simdata_list)

  return(list(ProteinLevelData=simdata))
}
