Converter_TPP<-function(x,CARRIER=FALSE){
  if(any(names(x)=="ProteinLevelData")){
    x<-x$ProteinLevelData
  }
  if(any(names(x)=="Channel")&any(names(x)=="temperature")){
    temps<-data.frame(temperature=as.character(unique(x$temperature)),Channel=as.character(unique(x$Channel)))
  }else if(any(names(x)=="temp_ref")){
    temps<-data.frame(temperature=as.character(unique(x$temperature)),Channel=as.character(unique(x$temp_ref)))
  }
  if(isTRUE(CARRIER)&any(names(x)=="Channel")){
    x$Channel<-ifelse(stringr::str_detect(x$Channel,"131N"),"131",x$Channel)
    x<-x|>dplyr::filter(Channel!="131C")
    }
  #Original data
  if(!any(names(x)=="sample_id")&any(names(x)=="Subject")){
    x$sample_id<-x$Subject
  }
  if(any(names(x)=="sample_id"&any(names(x)=="treatment"))){
    mappingMSstat<-x|>dplyr::select(Condition,sample_id,treatment,TechRepMixture)|>
      dplyr::distinct()|>
      dplyr::group_by(Condition,treatment)|>
      dplyr::mutate(Subject=paste0(Condition,"_",TechRepMixture),
                    Mixture=Subject,
                    Condition=as.character(Condition))
    x<-x|>dplyr::select(-Subject,-TechRepMixture,-treatment)|>
      dplyr::mutate(Condition=as.character(Condition))|>
                    dplyr::inner_join(mappingMSstat,by=c("Condition","sample_id"))

  }else if(any(stringr::str_detect(x$Mixture,"_[[:digit:]]+"))){

    x$TechRepMixture<-stringr::str_extract(stringr::str_extract(x$Mixture,"_[[:digit:]]+"),
                                           "[[:digit:]]+")
  }
  if(!any(names(x)=="vehicle")&length(unique(x$Condition))==2){
    x$Condition<-ifelse(stringr::str_detect(
      stringr::str_to_lower(x$treatment),
      "vehicle"),"Vehicle","Treatment")
    if(!any(names(x)=="Experiment")){
      x$treatment<-ifelse(stringr::str_to_lower(x$treatment)=="vehicle","Vehicle","Treatment")
      x$Experiment<-paste0(x$treatment,"_",x$TechRepMixture)
    }
    x$Mixture<-x$Experiment
  }else if(!any(names(x)=="Experiment")){
    x$Experiment<-paste0(x$treatment,"_",x$TechRepMixture)
    x$Mixture<-x$Experiment
  }
  if(any(isTRUE(class(x)))=="list"&any(stringr::str_detect(names(x),"TPPdata"))){
    x<-x$TPPdata
  }else if(any(names(x)=="ProteinLevelData")){
    if(stringr::str_detect(x$ProteinLevelData$Protein[1],"Sim")){
      x<-x$ProteinLevelData
      x$Mixture<-as.character(x$Experiment)
      x$Condition<-as.character(x$Condition)
    }else{
      x<-x$ProteinLevelData
    }
  }else if(!any(names(x)=="Mixture")&any(names(x)=="Experiment")){
    x<-x|>as.data.frame()
    x$Mixture<-x$Experiment
  }
  if(any(names(x)=="sample_id")){
    x$Subject<-x$sample_id
  }
  if(!any(names(x)=="Spectrum.File")|!any(names(x)=="Mixture")){
    x$Spectrum.File<-x$Mixture
  }
  if(nchar(as.character(x$TechRepMixture))[1]>1&all(stringr::str_ends(x$Mixture,"[[:digit:]]+_.raw"))|!any(names(x)=="TechRepMixture")){

    if(all(stringr::str_ends(x$Spectrum.File,"[[:digit:]]+_.raw"))){
      cols<-data.frame(Subject=unique(x$Subject),
                       Mixture=stringr::str_remove(unique(stringr::str_to_lower(x$Spectrum.File)),"[[:lower:]][[:digit:]]+_"))
      cols$TechRepMixture<-stringr::str_remove(stringr::str_extract(cols$Mixture,"[[:digit:]]+_.raw"),"_.raw")
      cols$Mixture<-stringr::str_extract(stringr::str_remove(cols$Mixture,"_.raw"),"[[:lower:]]+_[:digit:]")
    }else{
      cols<-data.frame(Subject=unique(x$Subject),
                       Mixture=stringr::str_extract(unique(stringr::str_to_lower(x$Spectrum.File)),"[[:lower:]]+[[:digit:]]+_"))

      cols$TechRepMixture<-stringr::str_extract(cols$Mixture,"[[:digit:]]+")
    }
    x<-x|>dplyr::select(-Mixture)|>dplyr::inner_join(cols)
  }

  if(length(unique(x$Mixture))==4&!any(names(x)=="Subject")){
    mix_x<-x|>dplyr::select(Mixture,Run)|>
      dplyr::distinct()|>
      dplyr::mutate(Subject=paste0("F",seq(1:dplyr::n())),
                    TechRepMixture=stringr::str_extract(Mixture,"[[:digit:]]+"))
    x<-x|>dplyr::inner_join(mix_x)
  }else if(!any(names(x)=="Subject")){
    mix_x<-x|>dplyr::select(Run)|>
      dplyr::distinct()|>
      dplyr::mutate(Subject=paste0("F",seq(1:dplyr::n())))
    x<-x|>dplyr::inner_join(mix_x)
  }



  stopifnot(any(stringr::str_detect(names(x),"Subject|Mixture")))

  if(!length(unique(x$Channel))==length(unique(temps$Channel))){
    stop("Inconsistency in temperature files.")
  }
  if(any(names(x)=="Channel")&any(names(x)=="temperature")){
    df.temps<-temps
    temps<-tidyr::pivot_wider(df.temps,names_from=Channel,values_from=temperature)
    x<-x|>dplyr::select(-temperature)|>dplyr::inner_join(df.temps)
  }else if(any(names(x)=="temp_ref")&any(names(x)=="temperature")){
    df.temps<-temps
    temps<-tidyr::pivot_wider(df.temps,names_from=temp_ref,values_from=temperature)
    x<-x|>dplyr::select(-temperature)|>dplyr::inner_join(df.temps)
  }else{
    df.temps<-temps
    temps<-tidyr::pivot_wider(df.temps,names_from=Channel,values_from=temperature)

  }

  if(isTRUE(CARRIER)&any(stringr::str_detect(x$Channel,"131N"))){
    x$Channel<-ifelse(stringr::str_detect(x$Channel,"131N"),"131",x$Channel)
    x<-x[!stringr::str_detect(x$Channel,"131C"),]

    temps<-temps[!stringr::str_detect(names(temps),"131C")]
    names(temps)<-stringr::str_replace(names(temps),"131N","131")
  }


  if(!any(names(x)=="Experiment")&any(names(x)=="treatment")){

    x1<-x|>
      dplyr::select(Subject,Mixture,treatment,TechRepMixture)|>
      dplyr::distinct()|>
      dplyr::arrange(treatment)
    if(any(stringr::str_detect(x$TechRepMixture,"4"))){
      reps<-x1|>
        dplyr::select(Subject,Mixture,treatment)|>
        dplyr::distinct()|>
        dplyr::arrange(treatment)|>
        dplyr::mutate(TechRepMixture=as.factor(rep(c(1,2),2)))
      x<-x1|>dplyr::select(-TechRepMixture)|>dplyr::inner_join(reps)
    }
    x$Experiment<-paste0(x$treatment,"_",x$TechRepMixture)
  }
  #rename_TPP
  if(any(names(x)=="Protein")&!any(names(x)=="uniqueID")){
    x$uniqueID<-x$Protein
  }
  if(any(names(x)=="uniqueID")&!any(names(x)=="gene_name")){
    x$gene_name<-x$uniqueID
  }
  if(any(names(x)=="Protein")&!any(names(x)=="gene_name")){
    x$gene_name<-x$Protein
  }
  if(any(names(x)=="Channel")&!any(names(x)=="temp_ref")){
    x$temp_ref<-x$Channel
  }
  if(any(is.na(x$Channel))&any(!is.na(x$temp_ref))){
    x$Channel<-x$temp_ref
  }
  if(any(names(x)=="value")&!any(names(x)=="I")){
    x$I<-x$value
  }
  if(any(names(x)=="Abundance")&!any(names(x)=="I")){
    x$I<-x$Abundance
  }

  if(any(names(x)=="Subject")&!any(names(x)=="sample_id")){
    x$sample_id<-x$Subject
  }

  if(any(names(x)=="Mixture")&!any(names(x)=="sample_id")){
    x$sample_id<-x$Mixture
  }

  if(any(names(x)=="Abundance")&!any(names(x)=="I")){
    x$I<-x$Abundance
  }else if(!any(names(x)=="Abundance")&!any(names(x)=="I")){
    x$I<-x$Intensity
  }
  if(any(names(x)=="Protein")){
    x$uniqueID<-x$Protein
    x$gene_name<-x$Protein
  }else if(any(names(x)=="ProteinName")){
    x$uniqueID<-x$ProteinName
    x$gene_name<-x$ProteinName
  }
  x$Subject<-x$Mixture

  return(x)
}
