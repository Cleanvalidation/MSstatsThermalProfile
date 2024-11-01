set_temps <- function(n,temperatures,protein_path,sample_mapping_name=NA){
  if(!is.logical(sample_mapping_name)){
    TMT<-readxl::read_xlsx(sample_mapping_name) %>%
      dplyr::rename("temp_ref"="TMT_label","temperature"="Temperature","sample_id"="Sample","sample_name"="MS_sample_number","time_point"="Time_point")
    TMT$time_point<-as.factor(TMT$time_point)

  }else{

    TMT<-data.frame(NA)
    if(n==10){
      TMT <- data.frame(Channel = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131'), temperature = temperatures, stringsAsFactors = FALSE)
    }else if (n==11){
      TMT <- data.frame(Channel = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131N','131C'), temperature = temperatures, stringsAsFactors = FALSE)
    }else if (n == 16){
      TMT <- data.frame(Channel = c('126', '127N', '127C', '128N', '128C', '129N','129C', '130N', '130C', '131N','131C','132N','132C','133N','133C','134N'), temperature = temperatures,stringsAsFactors = FALSE)
    }
  }
  return(TMT)
}
