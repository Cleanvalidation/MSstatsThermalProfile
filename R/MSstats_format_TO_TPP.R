#prepares MSstats groupComparisons output to TPP
#' @importFrom stringr str_extract
MSstats_format_TO_TPP<-function(summarisedProteins,temps,CARRIER=TRUE){
  if(!any(names(summarisedProteins)=="treatment")&!any(names(summarisedProteins)=="TechRepMixture")){
    summarisedProteins$treatment<-stringr::str_extract(summarisedProteins$Condition,'[[:lower:]]+')
    if(any(names(summarisedProteins)=="File.ID")){
    labels<-summarisedProteins|>dplyr::select(treatment,File.ID)|>
      dplyr::distinct()|>
      dplyr::group_by(treatment)|>
      dplyr::mutate(TechRepMixture=seq(1,dplyr::n()))
    }else{
      summarisedProteins<-summarisedProteins|>
        dplyr::group_by(Group)|>
        dplyr::mutate(TechRepMixture=as.factor(seq(dplyr::n())),
                      uniqueID=Protein)|>
        dplyr::ungroup()
      labels<-summarisedProteins|>dplyr::select(treatment,Subject,TechRepMixture)|>
        dplyr::distinct()|>
        dplyr::mutate(Mixture = as.factor(paste0(treatment,"_",TechRepMixture)))

    }
    if(!any(names(summarisedProteins)=="Experiment")){
      summarisedProteins<- summarisedProteins |> dplyr::right_join(labels) |>
        dplyr::mutate(Experiment = as.factor(paste0(treatment,"_",TechRepMixture)),
                      TechRepMixture = as.factor(TechRepMixture))
    }
  }

  if (any(names(summarisedProteins)== "Condition")){
    labels<-summarisedProteins|>dplyr::filter(!Condition=="Norm")
    labels<-labels|>
      dplyr::select(Condition,Mixture)|>
      dplyr::filter(!is.na(Condition))|>
      dplyr::distinct()|>
      dplyr::group_by(Condition)|>
      dplyr::mutate(TechRepMixture=as.factor(seq(1,dplyr::n())),
                    treatment=stringr::str_extract(as.character(Condition),"(?<=_)[[:lower:]]+|[[:lower:]]+"),
                    Condition = ifelse(treatment == "treated","treatment","vehicle"),
                    Experiment = paste0(Condition,"_",TechRepMixture))|>
      dplyr::distinct()
    summarisedProteins <- summarisedProteins |>
      dplyr::select(-Condition)|>
      dplyr::right_join(labels)|>
      dplyr::mutate(Condition = Experiment)
  }

  #rename protein columns into one format
  if(!any(names(summarisedProteins)=="uniqueID")&any(names(summarisedProteins)=="Accession")){
    summarisedProteins$uniqueID<-summarisedProteins$gene_name<-summarisedProteins$Accession
  }else if(!any(names(summarisedProteins)=="uniqueID")&any(names(summarisedProteins)=="Protein")){
    summarisedProteins$uniqueID<-summarisedProteins$gene_name<-summarisedProteins$Protein
  }

  #remove carrier channel assuming it is in 131C and rename 131N to 131
  if(any(stringr::str_detect(summarisedProteins$Channel,"131N"))&isTRUE(CARRIER)){
    summarisedProteins$Channel<-summarisedProteins$Channel
    summarisedProteins$Channel<-ifelse(stringr::str_detect(summarisedProteins$Channel,"131N"),"131",summarisedProteins$Channel)
    summarisedProteins<-summarisedProteins[!stringr::str_detect(summarisedProteins$Channel,"131C"),]
    if(any(class(summarisedProteins$Channel)=="factor")){
      summarisedProteins$Channel<-droplevels(summarisedProteins$Channel)
    }
  }else if(any(stringr::str_detect(summarisedProteins$Channel,"131C"))&isTRUE(CARRIER)){
    summarisedProteins$Channel<-summarisedProteins$Channel
    summarisedProteins<-summarisedProteins[!stringr::str_detect(summarisedProteins$Channel,"131C"),]
    summarisedProteins$Channel<-droplevels(summarisedProteins$Channel)
  }
  if(any(names(summarisedProteins)=="Channel")&!any(names(summarisedProteins)=="temperature")){
    summarisedProteins<-summarisedProteins|>dplyr::inner_join(temps)
  }

  x<-summarisedProteins

  x<-Converter_TPP(x,CARRIER=CARRIER)


  Orig_data<-dplyr::bind_rows(x) |>
    dplyr::select(uniqueID,gene_name,sample_id,TechRepMixture,I,Channel,Experiment)|>dplyr::distinct()


  Orig_data$Channel<-ifelse(stringr::str_detect(Orig_data$Channel,"C"),
                            stringr::str_replace(Orig_data$Channel,"C","H"),
                            Orig_data$Channel)
  Orig_data$Channel<-ifelse(stringr::str_detect(Orig_data$Channel,"N"),
                            stringr::str_replace(Orig_data$Channel,"N","L"),
                            Orig_data$Channel)
  Orig_data<-Orig_data |>
    dplyr::mutate(sample_id=as.character(sample_id))|>
    dplyr::filter(!is.na(Channel))|>
    dplyr::distinct()

  #pivot original data to wide and preface rel_fc to channels

  Data_pivot<-tryCatch(tidylog::pivot_wider(
    Orig_data,
    id_cols = c("uniqueID","gene_name","TechRepMixture","Experiment"),
    names_from = Channel,
    names_prefix = "rel_fc_",
    names_sep = "_",
    names_repair = "minimal",
    values_from = I,
    values_fn=mean

  ),error=function(e){
    tidyr::pivot_wider(
      Orig_data,
      #id_cols = c("uniqueID","gene_name","TechRepMixture","Experiment"),
      names_from = Channel,
      names_prefix = "rel_fc_",
      names_sep = "_",
      names_repair = "check_unique",
      values_from = I,
      values_fn=unique

    )
  })

  if(any(stringr::str_detect(names(Data_pivot),"rel_fc_NA"))){
    error("One of the channels is missing, please recheck the input")
  }
  pivot_x<-Data_pivot |> dplyr::distinct()
  check<-names(pivot_x)
  #select columns that contain channel info
  Ch_only<-check[stringr::str_detect(check,"[:digit:][:upper:]|[[:digit:]]+")]
  #column numbers that have reporter ion data
  cols_Ch_only<-which(check %in% Ch_only)


  #generate the configuration file from the original data
  config<-x |>
    dplyr::select(uniqueID,gene_name,Experiment,Condition,sample_id,TechRepMixture)|>
    dplyr::distinct()


  #join pivot data
  #config<-pivot_x |> dplyr::inner_join(config)
  config$treatment<-stringr::str_remove(stringr::str_extract(stringr::str_to_lower(config$Experiment),"[[:lower:]]+_"),"_")
  config$Condition<-ifelse(stringr::str_detect(stringr::str_to_lower(config$treatment),"vehicle"),"Vehicle","Treatment")

  config$ComparisonVT1<-NA
  config$ComparisonVT2<-NA

  config$ComparisonVT1<-ifelse(config$TechRepMixture==1,"x","")
  config$ComparisonVT2<-ifelse(config$TechRepMixture==2,"x","")

  df.temps<-data.frame(Channel=as.character(unique(x$Channel)),temperature=as.numeric(unique(x$temperature)))


  temps<-tidyr::pivot_wider(df.temps,names_from=Channel,values_from=temperature)
  temps<-furrr::future_imap(nrow(Data_pivot), ~temps)
  temps<-data.frame(temps[[1]])

  df.temps$Channel<-ifelse(stringr::str_detect(df.temps$Channel,"C"),
                           stringr::str_replace(df.temps$Channel,"C","H"),
                           df.temps$Channel)
  df.temps$Channel<-ifelse(stringr::str_detect(df.temps$Channel,"N"),
                           stringr::str_replace(df.temps$Channel,"N","L"),
                           df.temps$Channel)

  if(isTRUE(CARRIER)){
    df.temps$Channel<-ifelse(stringr::str_detect(df.temps$Channel,"131L"),"131",df.temps$Channel)
    df.temps<-df.temps[!stringr::str_detect(df.temps$Channel,"131H"),]
  }
  names(temps)<-df.temps$Channel
  Data_w_temps<-cbind(Data_pivot,temps)
  Data_w_temps<-Data_w_temps |> dplyr::inner_join(config)

  Data_w_temps$qssm<-NA
  Data_w_temps$qupm<-NA
  Data_w_temps$qssm<-as.integer(5)
  Data_w_temps$qupm<-as.integer(10)

  TPPconfig<-Data_w_temps |>
    dplyr::select(Experiment,ComparisonVT1,ComparisonVT2,dplyr::matches("[[:digit:]]+"),-dplyr::matches("rel_fc")) |>
    dplyr::distinct()|>
    dplyr::mutate(Experiment=as.character(Experiment)) |>
    dplyr::arrange(Experiment)

  TPPdata<-Data_w_temps |>
    dplyr::select(uniqueID,gene_name,Experiment,qssm,qupm,dplyr::matches("rel_fc_")) |>
    dplyr::distinct() |>
    dplyr::filter(!is.na(Experiment)) |>
    dplyr::arrange(Experiment) |>
    dplyr::mutate(Experiment=as.character(Experiment)) |>
    dplyr::group_split(gene_name)

  #TPPdata<-TPPdata |> purrr::keep(function(x) length(unique(x$Experiment))>=4)

  TPPdata<-dplyr::bind_rows(TPPdata) |> dplyr::group_split(Experiment)
  names(TPPdata)<-unique(TPPconfig$Experiment)

  resultPath<-file.path(getwd())

  TPPdata<-purrr::map(TPPdata,function(x) x |> dplyr::filter(!is.na(gene_name)))
  TPPdata<-lapply(TPPdata,function(x) x |> dplyr::mutate(uniqueID=as.factor(uniqueID),
                                                         gene_name=uniqueID))

  out<-list(TPPconfig,TPPdata)

  names(out)<-c("TPPconfig","TPPdata")
  return(out)
}

