AnnotationtoMSstatsTMT<-function(input,solvent="DMSO",temps){
  #pivot_longer
  input_long<-input %>%
    tidyr::pivot_longer(.,cols=names(input)[stringr::str_detect(names(input),"126|[[:digit:]]+N|[[:digit:]]+C")],
                        names_to="id",
                        values_to="value")

  input_long<-input_long %>%
    dplyr::mutate(temp_ref=stringr::str_extract(id,"126|[[:digit:]]+N|[[:digit:]]+C|[[:digit:]]+|131"))
  names(input_long)<-stringr::str_replace_all(names(input_long)," ","_")
  if(any(stringr::str_detect(input_long$Spectrum_File,"_[:digit:]"))){
    input_long$replicate<-stringr::str_extract(input_long$Spectrum_File,",_[:digit:]")
    input_long$replicate<-stringr::str_extract(input_long$Spectrum_File,"[:digit:]")
    input_long$TechRepMixture<-input_long$replicate
    input_long$replicate<-stringr::str_extract(input_long$Spectrum_File,",[:upper:]_[:digit:]")
    input_long$replicate<-stringr::str_extract(input_long$Spectrum_File,"[:digit:]")
    input_long$Mixture<-input_long$`File_ID`
  }else if(any(stringr::str_detect(input_long$Spectrum_File,".[:digit:]"))){
    input_long$replicate<-stringr::str_extract_all(input_long$id,".[:digit:]")
    input_long$replicate<-purrr::map(input_long$replicate,function(x)stringr::str_c(x,collapse=""))
    input_long$replicate<-purrr::map(input_long$replicate,function(x)as.factor(str_trunc(x,1,side=c("left"),ellipsis="")))
    input_long$replicate<-unlist(input_long$replicate)
    input_long<-input_long[!is.na(input_long$replicate),]
    input_long$TechRepMixture<-input_long$replicate
  }

  if(any(names(input_long)=="Intensity")){
    input_long<-input_long %>% dplyr::select(-Intensity)
  }

  if(any(names(input_long)=="Annotated_Sequence")){
    input_long<-input_long %>%
      dplyr::rename("ProteinName"="Master_Protein_Accessions",
                    "Run"="Spectrum_File",
                    #"PSM"="PSMs_Peptide_ID",
                    "Channel"="temp_ref",
                    "Abundance"="value") %>%
      dplyr::mutate(PrecursorCharge=Charge,
                    Mixture=Run,
                    BioReplicate=paste0(Run,"_",Channel),
                    Condition=paste0(Channel,"_",ifelse(stringr::str_detect(Run,solvent),"vehicle","treated")),
                    Group=Condition,
                    Experiment=File_ID)
  }else if(any(names(input_long)=="Master_Protein_Accessions")){
    input_long<-input_long %>%
      dplyr::rename("Protein"="Master_Protein_Accessions",
                    "Run"="Spectrum_File",
                    #"PSM"="PSMs_Peptide_ID",
                    "Channel"="temp_ref",
                    "Abundance"="value") %>%
      dplyr::mutate(PrecursorCharge=Charge,
                    Mixture=Run,
                    BioReplicate=paste0(Run,"_",Channel),
                    Condition=paste0(Channel,"_",ifelse(stringr::str_detect(Run,solvent),"vehicle","treated")),
                    Group=Condition,
                    Experiment=File_ID)
  }

  output<-input_long %>%
    dplyr::inner_join(temps) %>%
    dplyr::select(Protein,Run,Channel,Abundance,Mixture,BioReplicate)
  return(output)
}
