#add_variance_components
add_variance_components<-function(template,biovar,re_var,n_conditions,n_subjects,n_temps,t_range){
  modified_matrix <- matrix(NA, nrow = 20, ncol = n_conditions + 2)
  template<-template%>%dplyr::arrange(temperature)
  temperatures <- unique(template$temperature)[t_range]#get temperatures from the data
  for (i in 1:n_conditions) {
    for (j in 1:n_subjects) {
      sample_u <- rnorm(1, mean = 0, sd = sqrt(biovar))#simulate subject effect
      for (k in 1:n_temps) {
        current_responses <- template$Abundance[template$Condition == i][t_range[k]]
        sample_e <- rnorm(1, mean = 0, sd = sqrt(re_var))#simulate random error
        Y_true <- as.numeric(current_responses)
        Y_observed <- Y_true + sample_u + sample_e#add biological variance and error term
        #for a separate function, unlog Y_observed, average Y_obs for all temps, log
        modified_matrix[(j - 1)*n_temps + k , i] <- Y_observed
        modified_matrix[(j - 1)*n_temps + k , n_conditions + 1] <- temperatures[k]#set temperature
        modified_matrix[(j - 1)*n_temps + k , n_conditions + 2] <- j#set Subject
      }
    }
  }
  colnames(modified_matrix) <- c(seq(n_conditions), "temperature","Subject")
  modified_df <- as.data.frame(modified_matrix)
  #Prepare the annotation file data
  result_df_annotation<-prepare_annotation(template,modified_df)
  result_df_annotation$Simulation<-na.omit(result_df_annotation$Simulation)
  return(result_df_annotation$Simulation)
}
