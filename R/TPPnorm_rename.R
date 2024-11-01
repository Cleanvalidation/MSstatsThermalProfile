#normalizes data from the MSstats_TO_TPP output, disables FC filters
TPPnorm_rename<-function(input,df.temps=Channel2Temps_HumanData,filters=FALSE){

  tpptrInput<-input

  TRreqs<-TPP::tpptrDefaultNormReqs()
  if(!isTRUE(filters)){
  TRreqs$fcRequirements<-TRreqs$fcRequirements[1:2,]
  TRreqs$fcRequirements$thresholdLower<-c(-1,-1)
  TRreqs$fcRequirements$thresholdUpper<-c(Inf,Inf)
  }
  tpptrData<-TPP::tpptrNormalize(data=tpptrInput,normReqs = TRreqs)

  normTPP<-tpptrData$normData

  vault<-purrr::map(normTPP,function(x) Biobase::exprs(x))
  #use list names to generate experiment column, protein, techrepmixture and condition
  vault_names<-purrr::map2(vault,as.list(names(vault)),function(x,y)as.data.frame(x) %>%
                             dplyr::mutate(Experiment=y) %>%
                             dplyr::mutate(Protein=row.names(.),
                                           Accession=row.names(.),
                                           TechRepMixture=stringr::str_extract(Experiment,"[[:digit:]]+"),
                                           Condition=stringr::str_extract(stringr::str_to_lower(Experiment),"[[:lower:]]+")))
  #Bind data frames after generating columns
  vault_names<-dplyr::bind_rows(vault_names)
  vault_pivot<-vault_names%>%
    tidyr::pivot_longer(cols=names(vault_names)[stringr::str_detect(names(vault_names),"rel_fc")],
                        names_to="id",
                        values_to="Abundance")
  #convert channels back to C and N
  vault_rename<-vault_pivot %>%
    dplyr::mutate(Channel=stringr::str_extract(id,c('[:digit:][:digit:][:digit:][L|H]|126|131L|131H'))) %>%
    dplyr::mutate(Channel=stringr::str_replace(Channel,"H","C")) %>%
    dplyr::mutate(Channel=stringr::str_replace(Channel,"L","N")) %>%
    dplyr::select(-id) %>%
    dplyr::mutate(Condition=ifelse(Condition=="treatment","treated","vehicle"))
  col_info<-dplyr::bind_rows(vault_rename) %>%
    dplyr::select(Experiment,Condition,dplyr::starts_with("sample_id"),dplyr::starts_with("Spectrum")) %>%
    dplyr::mutate(Condition=ifelse(Condition=="Treatment","treated","vehicle")) %>%
    dplyr::distinct(.)
  vault_rename[] <- lapply( vault_rename, factor)
  vault_rename$Abundance<-as.numeric(as.character(vault_rename$Abundance))

  vault_join<-vault_rename %>% dplyr::left_join(col_info,by=c("Experiment","Condition"))
  if(any(names(vault_join)=="sample_id")){
  vault_join$Mixture<-vault_join$sample_id
  vault_join$BioReplicate<-vault_join$sample_id
  }else{#if this is a simulation, generate Mixture and bioReplicate
    df_names<-data.frame(Experiment=unique(vault_join$Experiment)) %>%
      dplyr::mutate(Mixture=paste0("F",rep(1:dplyr::n())),BioReplicate=Mixture)
    vault_join<-vault_join %>% dplyr::inner_join(df_names,by="Experiment")
  }
  if(any(names(vault_join)=="Spectrum.File")){
  vault_join$Run<-paste0(vault_join$Spectrum.File,"raw")
  }else{
    vault_join$Run<-"Simulation"
  }
  if(any(names(df.temps)=="temp_ref")){
  df.temps$Channel<-df.temps$temp_ref
  }else{
    df.temps$Channel<-as.factor(df.temps$Channel)
  }
  vault_join<-vault_join %>% dplyr::left_join(df.temps)
  vault_join$Condition<-paste0(vault_join$Channel,"_",vault_join$Condition)
  # Mixture, TechRepMixture, Run, Channel, Protein, Abundance, BioReplicate, Condition
  norm_TPP<-vault_join

  message(length(unique(norm_TPP$Protein))," proteins converted from TPP normalize to MSStatsTMT format.")

  return(norm_TPP)
}
