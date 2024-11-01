TPP_normalization<-function(pd_protein_data,TPPfilters,temps,reference,CARRIER=FALSE){
  if(any(names(pd_protein_data)=="Abundance")){
    pd_protein_data$value<-pd_protein_data$Abundance
  }
  #normalize data
  #if the condition is set to MSstats format channel_treatment
  if(any(stringr::str_detect(pd_protein_data$Condition,"[:punct:]"))){
    pd_protein_data<-pd_protein_data|>dplyr::group_by(Accession,File.ID,treatment)|>
      dplyr::mutate(value=suppressWarnings(value/mean(value[Channel==reference])),
                    uniqueID=Accession,
                    Subject=Experiment)|>
      dplyr::distinct()|>
      dplyr::ungroup()
  }else{ #if the condition is set to TPP format (e.g. vehicle or treated)
  pd_protein_data<-pd_protein_data|>dplyr::group_by(Accession,File.ID,Condition)|>
    dplyr::mutate(value=suppressWarnings(value/mean(value[Channel==reference])),
                  uniqueID=Accession,
                  Subject=Experiment)|>
    dplyr::distinct()|>
    dplyr::ungroup()
  }
  #rename data to MSStatsTMT format
  start=proc.time()
  Data<-MSstats_format_TO_TPP(pd_protein_data,temps=temps,CARRIER=CARRIER)
  end=proc.time()
  print(paste0("Renamed data to match TPP format in ",as.numeric(signif((end-start)[1],2))," seconds"))

  Data$TPPdata<-lapply(Data$TPPdata,function(x) x|>dplyr::select(-uniqueID)|>
                         dplyr::mutate(gene_name=as.character(gene_name)))
  #Check column order
  data(hdacTR_smallExample)
  order_names<-names(hdacTR_config)
  order_names<-stringr::str_replace(names(hdacTR_config),"131L","131")
  Data$TPPconfig<-Data$TPPconfig|>dplyr::mutate(Condition=stringr::str_extract(stringr::str_to_lower(Experiment),"[[:lower:]]+"))|>
                                                dplyr::select(order_names)
  #Replace vehicle and treatment with uppercase
  Data$TPPconfig<-Data$TPPconfig|>dplyr::mutate(Condition=ifelse(stringr::str_detect(Experiment,"vehicle"),"Vehicle","Treatment"),
                                                Experiment=ifelse(stringr::str_detect(Experiment,"vehicle"),
                                                                  stringr::str_replace(Experiment,"vehicle","Vehicle"),
                                                                  stringr::str_replace(Experiment,"treated","Treatment")))
  #Replace vehicle and treatment with uppercase
 Data$TPPdata<-lapply(Data$TPPdata,function(x) x|>
                        dplyr::mutate(Experiment=ifelse(stringr::str_detect(Experiment,"vehicle"),
                                                        stringr::str_replace(Experiment,"vehicle","Vehicle"),
                                                        stringr::str_replace(Experiment,"treated","Treatment"))))
 names(Data$TPPdata)<-unique(Data$TPPconfig$Experiment)
 Data$TPPconfig<-Data$TPPconfig|>dplyr::select(Experiment,Condition,ComparisonVT1,ComparisonVT2,"126","127L","127H","128L","128H","129L","129H","130L","130H","131")
  if(any(names(Data$TPPdata[[1]])=="uniqueID")){
    Data$TPPdata <- furrr::future_map(Data$TPPdata,function(x) x|>
                                        dplyr::select(-dplyr::starts_with("uniqueID"),
                                                      -dplyr::starts_with("variable"),
                                                      -dplyr::starts_with("value"))|>
                                        as.data.frame())
  }

  Data$TPPdata<-lapply(Data$TPPdata,function(x) tibble::column_to_rownames(x,"gene_name"))
  Data$TPPdata<-lapply(Data$TPPdata,function(x) x|>dplyr::mutate(gene_name=rownames(x)))

  order_names<-names(hdacTR_data$Vehicle_1)
  order_names<-stringr::str_replace(names(hdacTR_data$Vehicle_1),"131L","131")

  Data$TPPdata<-lapply(Data$TPPdata,function(x) x|>dplyr::select(order_names))

  tpptrData<-TPP::tpptrImport(configTable = Data$TPPconfig,data=Data$TPPdata,qualColName = "")
  TRreqs<-TPP::tpptrDefaultNormReqs()
  if(!isTRUE(TPPfilters)){

    TRreqs$fcRequirements<-TRreqs$fcRequirements
    TRreqs$fcRequirements$thresholdLower<-rep(0,3)
    TRreqs$fcRequirements$thresholdUpper<- rep(10,3)
    TRreqs$otherRequirements$thresholdLower<-1
  }


  normData<-TPP::tpptrNormalize(data=tpptrData,normReqs = TRreqs)

  ##################
  normData$normData<-TPP::tpptrTidyUpESets(normData$normData)
  normData$normData<- normData$normData |>dplyr::rename(temperature=x,
                                                        Abundance=y,
                                                        Channel=labelName,
                                                        Experiment=experiment)|>
    dplyr::mutate(Protein=uniqueID,gene_name=uniqueID,
                  Condition=stringr::str_extract(stringr::str_to_lower(Experiment),
                                                 "[[:lower:]]+"))|>
    dplyr::mutate(Channel=stringr::str_replace(Channel,"H","C")) |>
    dplyr::mutate(Channel=stringr::str_replace(Channel,"L","N"),
                  Condition=ifelse(stringr::str_detect(Condition,"treat"),"treated","vehicle"))

  normData$normData<-normData$normData|>dplyr::mutate(Condition=paste0(Channel,"_",Condition),
                                                      TechRepMixture=stringr::str_extract(replicate,"[[:digit:]]+"))|>
    dplyr::mutate(Condition=ifelse(stringr::str_detect(Channel,reference),"Norm",Condition),
                  Mixture=Experiment,
                  BioReplicate=ifelse(!any(stringr::str_detect(names(normData$normData),"Run|Spectrum.File")),1,NA))
  if(!any(names(normData$normData)=="Run")){
    normData$normData$Run<-normData$normData$Mixture
  }
  all_x<-data.frame(Experiment=unique(normData$normData$Experiment))|>dplyr::mutate(Subject=unique(normData$normData$Mixture))
  normData$normData<-normData$normData|>dplyr::inner_join(all_x)
  normData$normData$BioReplicate<-normData$normData$Mixture
  return(normData)
}

