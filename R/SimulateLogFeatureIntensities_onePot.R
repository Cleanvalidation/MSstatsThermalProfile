SimulateLogFeatureIntensities_onePot<-function(Yobs_crt)
{

  #Selecting the required columns from simulation result
  template_protein<-data.frame(Protein=Yobs_crt$Protein)
  template_protein$temperature<-Yobs_crt$temperature
  template_protein$Abundance<-Yobs_crt$Abundance
  template_protein$cond_repl<-paste0(Yobs_crt$treatment,"_",stringr::str_split_fixed(Yobs_crt$Subject,"_",2)[,2])
  # Making a matrix of abundances for each condition and replicate
  data<-template_protein %>%
    group_by(cond_repl ) %>%
    mutate(row = row_number()) %>%
    tidyr::pivot_wider(names_from = cond_repl, values_from = Abundance) %>%
    select(-row)
  colnames<- unique(template_protein$cond_repl)
  protein_matrix <- as.matrix(data[, colnames])
  # Adapt column and row maxquant_protein_table
  colnames(protein_matrix) <- colnames
  rownames(protein_matrix) <- data$Protein
  # Unloging, averaging, loging back for one pot scenario
  matrix<-protein_matrix
  #Unlogging
  combined_matrix<-2^(matrix)
  #Averaging
  combined_matrix <- as.matrix(aggregate(combined_matrix, by = list(rownames(matrix)), FUN = mean,na.rm=TRUE))
  colnames <- grep("\\_", colnames(combined_matrix), value=TRUE)
  #Generating a dataframe from matrix
  df<-as.data.frame(combined_matrix)
  df<-df %>% pivot_longer(cols=colnames,
                          names_to='cond_rep',
                          values_to='Abundance')
  #loging back
  df$Abundance<-log2(as.numeric(df$Abundance))
  #Adding required columns for msstatscomparison
  df$Condition<-stringr::str_split_fixed(df$cond_rep,"_",2)[,1]
  df$BioReplicate<-stringr::str_split_fixed(df$cond_rep,"_",2)[,2]
  colnames(df)[1] <- "Protein"
  df$Run<-df$cond_rep
  df$Channel<-rep("128",nrow(df))
  df$TechRepMixture<-rep(1,nrow(df))
  df$Mixture<-df$cond_rep
  return(df)
}
