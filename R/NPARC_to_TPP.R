NPARC_to_TPP<-function(x,temps=df.temps,string=FALSE,NORM=FALSE,TPPnormfilters=FALSE){#rename script data to run TPP
  if(any(names(x)=="Accession")&!any(names(x)=="uniqueID")){
    x$uniqueID<-x$Accession
  }
  if(any(names(x)=="uniqueID")&!any(names(x)=="gene_name")){
    x$gene_name<-x$uniqueID
  }
  if(any(names(x)=="treatment")&!any(names(x)=="Condition")){
    x$Condition<-x$treatment
  }
  if(any(names(x)=="dataset")&!any(names(x)=="Condition")){
    x$Condition<-x$dataset
  }
  if(any(names(x)=="uniqueID")&!any(names(x)=="gene_name")){
    x$gene_name<-x$uniqueID
  }
  if(any(names(x)=="C")&!any(names(x)=="temperature")){
    x$temperature<-x$C
  }
  if(any(names(x)=="value")&!any(names(x)=="I")){
    x$I<-x$value
  }
  if(any(names(x)=="relAbundance")&!any(names(x)=="I")){
    x$I<-x$relAbundance
  }
  #add replicate columns
  if(!any(names(x)=="replicate")){
    x<-replicate_labels(x)
  }
  data<-dplyr::bind_rows(x)
  if(!any(names(data)=="Spectrum.File")&any(names(data)=="compoundConcentration")){
    data$dataset<-data$Condition<-ifelse(data$compoundConcentration==0,"Vehicle","Treatment")
    data<-data|>dplyr::mutate(sample_id=paste0(dataset,"_",replicate))
    data$Spectrum.File<-paste0(data$sample_id,".raw")
  }else{
    data$replicate<-stringr::str_extract(data$Spectrum.File,"[[:upper:]][[:digit:]]+_|[[:lower:]][[:digit:]]+_")
    data$replicate<-stringr::str_extract(data$replicate,"[[:digit:]]+")
  }

  data<-data|>dplyr::inner_join(temps)

  data$temp_ref<-data$Channel

  data<-data |>
    dplyr::select(gene_name,Spectrum.File,sample_id,replicate,Condition,I,temp_ref) |>
    dplyr::filter(!is.na(I))

  data_wide<-pivot_wider(
    data,
    id_cols = NULL,
    names_from = temp_ref,
    names_prefix = "rel_fc_",
    names_sep = "_",
    names_repair = "minimal",
    values_from = I,
    values_fn=unique

  )

  data_wide<-data_wide %>% distinct(.)
  check<-names(data_wide)
  check1<-check[stringr::str_detect(check,"[:digit:][:upper:]|[[:digit:]]+")]
  #column numbers that have reporter ion data
  data2<-which(check %in% check1)
  #replace C or N with L and H
  check1<-stringr::str_replace(check1,"C","H")
  check1<-stringr::str_replace(check1,"N","L")
  #replace names
  check[data2]<-check1
  names(data_wide)<-check
  data_wide$Condition<-ifelse(data_wide$Condition=="Vehicle","Vehicle","Treatment")
  if(any(stringr::str_detect(stringr::str_to_lower(data_wide$sample_id),"[[:lower:]]+_[[:digit:]]+"))){
  data_wide$Experiment<-data_wide$sample_id
  }else{
  data_wide$Experiment<-ifelse(data_wide$Condition=="Treatment",
                               paste0("Treatment_",data_wide$replicate),paste0(data_wide$Condition,"_",data_wide$replicate))
  }
  data_wide$ComparisonVT1<-NA
  data_wide$ComparisonVT2<-NA

  data_wide$ComparisonVT1<-ifelse(data_wide$replicate==1,"x","")
  data_wide$ComparisonVT2<-ifelse(data_wide$replicate==2,"x","")

  check1<-check[stringr::str_detect(check,"rel_fc_[[:digit:]]+|rel_fc_[[:digit:]]+[:upper:]")]
  #column numbers that have reporter ion data
  data2<-which(check %in% check1)

  config<-dplyr::bind_rows(data_wide)

  check<-c(config %>% dplyr::select(Experiment,Condition,ComparisonVT1,ComparisonVT2),config[data2])
  temp_ref<-stringr::str_replace(temps$temp_ref,"C","H")
  temp_ref<-stringr::str_replace(temp_ref,"N","L")


  temps<-tidyr::pivot_wider(temps,names_from=Channel,values_from=temperature,values_fn = "unique")
  temps<-purrr::map_dfr(seq_len(nrow(data)), ~temps)
  temp_ref<-stringr::str_replace(names(temps),"C","H")
  temp_ref<-stringr::str_replace(temp_ref,"N","L")

  names(temps)<-temp_ref
  data_wide<-cbind(data_wide,temps[1,])

  # data_wide$qssm<-as.integer(sample_id(0:125,nrow(data_wide),replace=TRUE))
  # data_wide$qupm<-as.integer(sample_id(4:40,nrow(data_wide),replace=TRUE))
  data_wide$qssm<-as.integer(5)
  data_wide$qupm<-as.integer(10)


  config<-dplyr::bind_rows(data_wide)|>
    dplyr::select(Experiment,Condition,ComparisonVT1,ComparisonVT2,dplyr::starts_with("1"))|>
    dplyr::distinct()
  #data_wide<-data_wide |> dplyr::mutate(gene_name = gsub("_IPI.*", "", gene_name))
  data_wide$gene_name<-stringr::str_replace_all(data_wide$gene_name," ",".")
  #data_wide$gene_name<-make.unique(data_wide$gene_name)
  data_wide$gene_name<-as.factor(data_wide$gene_name)

  if(isTRUE(NORM)){
    data_wide<-TPPnorm_rename(data_wide,df.temps=temps,filters=TPPnormfilters)
    data_wide<-data_wide|>
      dplyr::select(gene_name,Experiment,sample_id,qssm,qupm,dplyr::starts_with("rel"))

    data_wide<-data_wide|>dplyr::group_by(sample_id)|>dplyr::group_split()|>lapply(function(x) x[,order(names(x))])

    data_wide<-lapply(data_wide,function(x) x|>dplyr::select(-sample_id)|>distinct())
    names(data_wide)<-unique(config$Experiment)
    data_wide=list(TPPconfig=config,TPPdata=data_wide)
    return(data_wide)
  }
  data_wide<-data_wide|>
    dplyr::select(gene_name,Experiment,sample_id,qssm,qupm,dplyr::starts_with("rel"))

  rownames(data_wide)<-make.unique(as.character(data_wide$gene_name))

  data_wide<-data_wide|>dplyr::group_by(sample_id)|>dplyr::group_split()|>lapply(function(x) x[,order(names(x))])
  data_wide<-lapply(data_wide,function(x) x|>dplyr::select(-sample_id)|>distinct())
  names(data_wide)<-unique(config$Experiment)

  data_wide=list(TPPconfig=config,TPPdata=data_wide)
  return(data_wide)
}
