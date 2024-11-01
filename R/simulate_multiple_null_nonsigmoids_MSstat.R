
simulate_multiple_null_nonsigmoids_MSstat = function(runs=1000,n=4,error_sd=0.05,subject_sd=0,half_temp_grid=FALSE){
  simdata_list = vector(mode="list",length=runs)


  for(i in 1:runs){
    Protein = paste("Null_SubSD",subject_sd,"Sim",i,sep="_")

    sim = simulate_null_nonsigmoid_MSstat(n=n,error_sd=error_sd,subject_sd = subject_sd,half_temp_grid=half_temp_grid)
    simdata_list[[i]] = sim
    simdata_list[[i]]$Protein = Protein

  }

  simdata = dplyr::bind_rows(simdata_list)


  return(simdata)
}
