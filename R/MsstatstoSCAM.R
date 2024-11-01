#convert MSstats output to fit scam models
MSstatstoSCAM<-function(Result){
  if(any(names(Result)=="Condition")&!any(names(Result)=="treatment")){
    Result$treatment<-as.factor(Result$Condition)
  }
  if(any(names(Result)=="Protein")&!any(names(Result)=="Accession")){
    Result$Accession<-Result$Protein
  }
  if(!any(names(Result)=="Subject")){
    cols_subject<-data.frame(Mixture=unique(Result$Mixture))|>dplyr::mutate(Subject=paste0("F",seq(length(unique(Result$Mixture)))))
    Result<-Result|>dplyr::inner_join(cols_subject)
     }
  if(any(stringr::str_detect(Result$Condition[1],"_"))& !any(names(Result)=="treatment")){
    if(any(stringr::str_detect(Result$treatment[1],"_"))){
      Result<-Result |> tidyr::separate(Condition,into=c("temperature","treatment"),sep="_") |> dplyr::mutate(I=Abundance,
                                                                                                              replicate=TechRepMixture,
                                                                                                              Accession=Protein,
                                                                                                              Spectrum.File=Run,
                                                                                                              sample_id=Subject,
                                                                                                              Subject=Subject,
                                                                                                              Condition=treatment)|>
        dplyr::mutate(treatment=as.factor(treatment),uniqueID=Protein)
    }else{
  Result<-Result |> dplyr::mutate(I=Abundance,
                                                                                                  replicate=TechRepMixture,
                                                                                                  Accession=Protein,
                                                                                                  Spectrum.File=Run,
                                                                                                  sample_id=Subject,
                                                                                                  Subject=Subject,
                                                                                                  Condition=treatment)|>
    dplyr::mutate(treatment=as.factor(treatment),uniqueID=Protein)
    }
  }else{
    Result<-Result |>  dplyr::mutate(I=Abundance,
                                                                                                              replicate=TechRepMixture,
                                                                                                              Accession=Protein,
                                                                                                              Spectrum.File=Run,
                                                                                                              sample_id=Subject,
                                                                                                              Subject=Subject)|>
      dplyr::mutate(Condition=treatment,uniqueID=Protein)
  }
  if(any(names(Result)=="Experiment")){
  if(stringr::str_detect(Result$Experiment[1],"_")){
    Result$Condition<-Result$treatment<-as.factor(stringr::str_extract(stringr::str_to_lower(Result$Experiment),"[[:lower:]]+"))
  }
  }else{
    if(any(is.na(Result$Condition))&any(!is.na(Result$treatment))){
      Result$Condition<-as.factor(stringr::str_extract(stringr::str_to_lower(Result$treatment),"[[:lower:]]+"))
    }

  }
  if(any(stringr::str_detect(Result$Condition,"orm"))){
    Result$Condition<-Result$treatment<-as.factor(ifelse(stringr::str_extract(stringr::str_to_lower(Result$Experiment),"[[:lower:]]+")=="treatment","treated","vehicle"))
  }
  return(Result)
}
