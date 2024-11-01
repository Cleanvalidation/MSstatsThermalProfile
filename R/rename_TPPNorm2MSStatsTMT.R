#normTPP is the output from the tpptrNormalize function.  This function converts the Bioexprs objects to data frames
rename_TPPNorm2MSStatsTMT<-function(normTPP){
  vault<-purrr::map(normTPP,function(x) Biobase::exprs(x))
  vault_names<-purrr::map2(vault,as.list(names(vault)),function(x,y)as.data.frame(x) %>%
                             dplyr::mutate(Experiment=y) %>%
                             dplyr::mutate(Protein=row.names(.),
                                           Accession=row.names(.),
                                           TechRepMixture=stringr::str_extract(Experiment,"[[:digit:]]+"),
                                           Condition=stringr::str_extract(stringr::str_to_lower(Experiment),"[[:lower:]]+")))
  vault_names<-dplyr::bind_rows(vault_names)
  vault_pivot<-vault_names%>%
    tidyr::pivot_longer(cols=names(vault_names)[stringr::str_detect(names(vault_names),"rel_fc")],
                        names_to="id",
                        values_to="Abundance")
  vault_rename<-vault_pivot %>%
    dplyr::mutate(Channel=stringr::str_extract(id,c('[:digit:][:digit:][:digit:][L|H]|126|131L|131H'))) %>%
    dplyr::mutate(Channel=stringr::str_replace(Channel,"H","C")) %>%
    dplyr::mutate(Channel=stringr::str_replace(Channel,"L","N")) %>%
    dplyr::select(-id) %>%
    dplyr::mutate(Condition=ifelse(Condition=="treatment","treated","vehicle"))
  col_info<-dplyr::bind_rows(TPP_Cliff) %>%
    tidyr::pivot_longer(cols=names(TPP_Cliff)[stringr::str_detect(names(TPP_Cliff),"rel_fc")],
                        names_to="id",
                        values_to="Abundance") %>%
    dplyr::select(Experiment,Condition,sample_id,dplyr::starts_with("Spectrum")) %>%
    dplyr::mutate(Condition=ifelse(Condition=="Treatment","treated","vehicle")) %>%
    distinct(.)
  vault_rename[] <- lapply( vault_rename, factor)
  vault_rename$Abundance<-as.numeric(as.character(vault_rename$Abundance))

  vault_join<-vault_rename %>% dplyr::left_join(col_info,by=c("Experiment","Condition"))
  vault_join$Mixture<-vault_join$sample_id
  vault_join$BioReplicate<-vault_join$sample_id

  if(any(names(vault_join)=="Spectrum.File")){
    vault_join$Run<-paste0(vault_join$Spectrum.File,"raw")
  }else{
    vault_join$Run<-"Simulation"
  }
  df.temps$Channel<-df.temps$temp_ref
  vault_join<-vault_join %>% dplyr::left_join(df.temps)
  vault_join$Condition<-paste0(vault_join$temp_ref,"_",vault_join$Condition)
  # Mixture, TechRepMixture, Run, Channel, Protein, Abundance, BioReplicate, Condition
  norm_MSStats_format<-vault_join

  message(length(unique(norm_MSStats_format$Protein))," proteins converted from TPP normalize to MSStatsTMT format.")
  return(norm_MSStats_format)
}
