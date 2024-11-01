replicate_labels<-function(x){
  if(any(names(x)=="Spectrum.File")&any(stringr::str_detect(x$Spectrum.File,"[:upper:][:digit:]+|[:lower:][[:digit:]]+"))){
    if(any(stringr::str_detect(x$Spectrum.File,"[:upper:][:digit:]+|[:lower:][[:digit:]]+"))){
      x$replicate<-stringr::str_extract(x$Spectrum.File,"[:upper:][:digit:]+|[:lower:][[:digit:]]+")
      x$replicate<-as.factor(stringr::str_extract(x$replicate,"[[:digit:]]+"))
      return(x)
    }
  }else if(any(names(x)=="Spectrum_File")&any(stringr::str_detect(x$Spectrum_File,"[:upper:][:digit:]+|[:lower:][[:digit:]]+"))){
    if(any(stringr::str_detect(x$Spectrum_File,"[:upper:][:digit:]+|[:lower:][[:digit:]]+"))){
      x$replicate<-stringr::str_extract(x$Spectrum_File,"[:upper:][:digit:]+|[:lower:][[:digit:]]+")
      x$replicate<-as.factor(stringr::str_extract(x$replicate,"[[:digit:]]+"))
      return(x)
    }
  }
  if(any(stringr::str_detect(names(x)," "))){
    names(x)<-stringr::str_replace_all(names(x)," ","_")
  }
  #add replicates
  if(any(names(x)=="Annotated_Sequence")&!(any(names(x)=="replicate"))){
    if(any(names(x)=="treatment")){
      x<-dplyr::bind_rows(x) %>%
        distinct(.) %>%
        dplyr::group_by(uniqueID,Annotated_Sequence,temp_ref,treatment,sample_name) %>%
        dplyr::group_split()
      x<-purrr::map(x,function(x)x %>% dplyr::mutate(Replicate=row.names(.),
                                                     replicate=row.names(.)))
    }else{
      x<-dplyr::bind_rows(x) %>%
        distinct(.) %>%
        dplyr::mutate(replicate=
                        stringr::str_extract(Spectrum_File,"[[:digit:]]+."))

      x<-x %>% dplyr::mutate(replicate=stringr::str_extract(replicate,"[[:digit:]]+"))

    }
  }else if (any(names(x)=="Accession")&any(names(x)=="replicate")){
    x<-dplyr::bind_rows(x)
  }else if(!any(names(x)=="replicate")&any(names(x)=="uniqueID")){
    x<-dplyr::bind_rows(x) %>%
      distinct(.) %>%
      dplyr::group_by(uniqueID,temp_ref,treatment,sample_name) %>%
      dplyr::group_split()
    x<-purrr::map(x,function(x)x %>% dplyr::mutate(Replicate=row.names(.),
                                                   replicate=row.names(.)))
  }else if(any(names(x)=="temperature")&any(names(x)=="Accession")){
    x<-dplyr::bind_rows(x) %>%
      distinct(.) %>%
      dplyr::group_by(Accession,temperature,treatment,sample_name) %>%
      dplyr::group_split()
    x<-purrr::map(x,function(x)x %>% dplyr::mutate(Replicate=row.names(.),
                                                   replicate=row.names(.)))
  }else if(any(names(x)=="Accession")&any(stringr::str_detect(names(x),"Found"))){
    x<-x[,!stringr::str_detect(names(x),"Found")]

  }else if(any(names(x)=="Mixture")&!any(names(x)=="temp_ref")){
    x<-dplyr::bind_rows(x) %>%
      distinct(.) %>%
      dplyr::group_by(Accession,temp_ref,Mixture) %>%
      dplyr::group_split()
    x<-purrr::map(x,function(x)x %>% dplyr::mutate(Replicate=row.names(.),
                                                   replicate=row.names(.)))

  }
}
