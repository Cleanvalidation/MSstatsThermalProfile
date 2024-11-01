MSstatsTMT_Proteinconverter<-function(df_raw,solvent="DMSO",ref="126",Frac=TRUE,CFS=FALSE,CARRIER=FALSE){
  editPSMs2<-data.frame()
  ref<-as.character(ref)

  input_files<-getwd()

  #if there are no technical replicates, add them
  editPSMs2<-replicate_labels(df_raw)
  editPSMs2<-dplyr::bind_rows(editPSMs2)

  if(!any(names(editPSMs2)=="replicate")){
    editPSMs2<-replicate_labels(editPSMs2)
  }
  if(any(names(editPSMs2)=="Intensity")){
    editPSMs2<-editPSMs2 %>% dplyr::select(-Intensity)
  }

  #PeptideSequence, PrecursorCharge, FragmentIon and ProductCharge
  if(any(names(editPSMs2)=="Spectrum.File")&any(names(editPSMs2)=="I")){
    editPSMs2$Run<-paste0(ifelse(stringr::str_detect(editPSMs2$Spectrum.File,"NOcarrier")==TRUE,"nC",ifelse(stringr::str_detect(editPSMs2$Spectrum.File,"carrier")==TRUE,"C",NA)),'_',
                          ifelse(stringr::str_detect(editPSMs2$Spectrum.File,"NO_FAIMS")==TRUE,"nF",ifelse(stringr::str_detect(editPSMs2$Spectrum.File,"r_FAIMS")==TRUE,"F",NA)),'_',
                          ifelse(stringr::str_detect(editPSMs2$Spectrum.File,"S_eFT")==TRUE,"E",ifelse(stringr::str_detect(editPSMs2$Spectrum.File,"S_Phi")==TRUE,"S",NA)))
    editPSMs2<-editPSMs2 %>%
      dplyr::rename("Protein"="Accession",
                    "Abundance"="I",
                    "Channel"="temp_ref",
                    "Subject"="sample_id") %>%
      dplyr::mutate(PeptideSequence=NA,
                    PrecursorCharge=NA,
                    FragmentIon=NA,
                    ProductCharge=NA,
                    BioReplicate=paste0(Run,"_",Channel),
                    Mixture=ifelse(stringr::str_detect(Spectrum.File,solvent),
                                   paste0(Channel,"_","vehicle"),
                                   paste0(Channel,"_","treated")),
                    Group=Mixture,
                    Condition=ifelse(stringr::str_detect(Spectrum.File,solvent),
                                     paste0(Channel,"_","vehicle"),
                                     paste0(Channel,"_","treated")),
                    TechRepMixture=replicate)
    if(all(is.na(editPSMs2$Spectrum.File))){
      editPSMs2$Condition=NA
    }else{
      editPSMs2$Condition<-as.factor(editPSMs2$Condition)
    }
  }

   df_raw<-dplyr::bind_rows(df_raw)

  annotation<-editPSMs2 %>%
    dplyr::select(Run,Fraction,TechRepMixture,Mixture,Channel,BioReplicate,Condition) %>%
    distinct(.)
  if(any(stringr::str_detect(editPSMs2$Channel,ref))){
    editPSMs2$Condition<-ifelse(stringr::str_detect(editPSMs2$Channel,ref),"Norm",as.character(editPSMs2$Condition))
  }
  # if(isTRUE(CARRIER)&any(stringr::str_detect(editPSMs2$Condition,"131N"))){
  #   editPSMs2$Condition<-ifelse(stringr::str_detect(editPSMs2$Condition,"131N"),stringr::str_replace(editPSMs2$Condition,"N",""),editPSMs2$Condition)
  # }
  if(any(!isFALSE(names(annotation) %in% c("Run","Fraction","TechRepMixture","Mixture","Channel","BioReplicate","Condition")))){
    return(list(ProteinLevelData=editPSMs2))
  }else{
    warning(paste0("Columns ",names(annotation)[isFALSE(names(annotation) %in% c("Run","Fraction","TechRepMixture","Mixture","Channel","BioReplicate","Condition"))]," are missing"))
  }

}
