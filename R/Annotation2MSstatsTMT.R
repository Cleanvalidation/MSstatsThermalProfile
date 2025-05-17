#' @importFrom stringr str_extract

Annotation2MSstatsTMT<-function(input,solvent="DMSO",temps,reference,CARRIER=TRUE){


  if(any(stringr::str_detect(names(input)," "))){
    names(input)<-stringr::str_replace_all(names(input)," ",".")
  }
  if(isTRUE(CARRIER)){
    if(any(stringr::str_detect(names(input),"131N"))){
      names(input)<-stringr::str_replace_all(names(input),"131N","131")
      input<-input[,!stringr::str_detect(names(input),"131C")]
    }
  }

  #pivot_longer
  input_long<-input |>
    tidyr::pivot_longer(cols=names(input)[stringr::str_detect(names(input),"[:digit:][:digit:][:digit:]|[[:digit:]]+N|[[:digit:]]+C")],
                        names_to="id",
                        values_to="Abundance")

  input_long<-input_long |>
    dplyr::mutate(Channel=stringr::str_extract(id,"126|[[:digit:]]+N|[[:digit:]]+C|131"))


  input_long$Subject<-stringr::str_extract(input_long$Spectrum.File,"_[[:upper:]]+[:digit:]_|[[:upper:]][[:lower:]]+[[:digit:]]+_")
  input_long$Mixture<-input_long$Subject
  input_long$replicate<-as.character(stringr::str_extract(input_long$Subject,"[[:digit:]]+"))
  input_long$TechRepMixture<-input_long$replicate

  if(nchar(input_long$TechRepMixture[1])>2|any(is.na(input_long$TechRepMixture))){#spectrum file columns can have different formats for subject, mixture and bioreplicate
    input_long$Subject<-stringr::str_extract(input_long$Spectrum.File,"_[[:upper:]]+[[:digit:]]+_|[[:upper:]][[:lower:]]+[[:digit:]]+_")
    input_long$Mixture<-input_long$Subject
    input_long$TechRepMixture<-as.character(stringr::str_extract(input_long$Mixture,"[[:digit:]]+"))
    input_long$replicate<-input_long$TechRepMixture
    if(nchar(input_long$TechRepMixture[1])>2|any(is.na(input_long$TechRepMixture))){
      input_long$Subject<-stringr::str_extract(stringr::str_to_lower(input_long$Spectrum.File),"[[:lower:]]+[[:digit:]]+")
      input_long$Mixture<-input_long$Subject
      input_long$TechRepMixture<-as.character(stringr::str_extract(input_long$Mixture,"[[:digit:]]+"))
      input_long$replicate<-input_long$TechRepMixture
      if(nchar(input_long$TechRepMixture[1])>2|any(is.na(input_long$TechRepMixture))){
        input_long$Subject<-stringr::str_extract(input_long$Spectrum.File,"[[:upper:]][[:lower:]]+_[[:digit:]]+")
        input_long$Mixture<-input_long$Subject
        input_long$TechRepMixture<-as.character(stringr::str_extract(input_long$Mixture,"[[:digit:]]+"))
        input_long$replicate<-input_long$TechRepMixture

      }
    }

  }
  #specific to Proteome Discoverer to extract fractions
  if(any(stringr::str_detect(input_long$File.ID,"."))){
    input_long$Fraction<-stringr::str_remove(input_long$File.ID,"F[[:digit:]]+.")
  }

    if(any(names(input_long)=="Intensity")){
    input_long<-input_long |> dplyr::select(-Intensity)
  }

  if(!is.na(reference)){
    if(any(stringr::str_detect(names(input_long),"Master"))&!any(stringr::str_detect(names(input_long),"Condition"))){
      #if only the master protein accessions is available
      input_long<-input_long |>
        dplyr::rename("Protein"="Master.Protein.Accessions",
                      "Run"="Spectrum.File") |>
        dplyr::mutate(Condition=ifelse(stringr::str_detect(Channel,reference),"Norm",paste0(Channel,"_",ifelse(stringr::str_detect(Run,solvent),"vehicle","treated"))),
                      BioReplicate=Mixture,
                      treatment=paste0(ifelse(stringr::str_detect(Run,solvent),"vehicle","treated")))
    }else if (any(names(input_long)=="Accession")){
      input_long<-input_long |>
        dplyr::rename("Protein"="Accession",
                      "Run"="Spectrum.File") |>
        dplyr::mutate(Condition=ifelse(stringr::str_detect(Channel,reference),"Norm",paste0(Channel,"_",ifelse(stringr::str_detect(Run,solvent),"vehicle","treated"))),
                      BioReplicate=Mixture,
                      treatment=paste0(ifelse(stringr::str_detect(Run,solvent),"vehicle","treated")))
    }

  }
  input_long$TechRepMixture<-as.character(input_long$TechRepMixture)
  if(nchar(input_long$TechRepMixture[1])>1){
    reps<-data.frame(Subject=as.character(unique(input_long$Subject)),TechRepMixture=as.character(seq(1:length(unique(input_long$Subject)))))
    input_long<-input_long|>dplyr::select(-TechRepMixture)|>dplyr::inner_join(reps)

  }
  class(input_long)<-"data.frame"
  if(any(stringr::str_detect(names(input_long),"Fraction"))){#if this study has fractions
    output<-input_long |>
      dplyr::inner_join(temps) |>
      dplyr::select(Run,Fraction,TechRepMixture,Channel,Subject,Condition,Mixture,BioReplicate,treatment) |>
      dplyr::distinct()
  }else{
    output<-input_long |>
      dplyr::inner_join(temps) |>
      dplyr::select(Run,TechRepMixture,Channel,Subject,Condition,Mixture,BioReplicate,treatment) |>
      dplyr::distinct()
  }
  stopifnot(nrow(output)>0)
  if(any(stringr::str_detect(output$Fraction,"F"))){
    output<-output|>dplyr::select(-Fraction)
  }
  if(!any(stringr::str_detect(output$Fraction,output$Run))){
    output$Fraction<-stringr::str_extract(stringr::str_extract(output$Run,"[[:digit:]]+.raw"),"[[:digit:]]+")
  }

  if(isTRUE(CARRIER)){
    return(list(df=input,annotation=output))
  }
  #output$BioReplicate<-paste0(output$Mixture,"_",output$Condition)
  return(output)
}

