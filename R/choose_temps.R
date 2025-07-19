choose_temps<-function(df,temps=NA){
  if(any(is.na(temps))){
    if(any(df$Run=="OnePot")){
      return(df)
    }else{
      temps<-unique(df$temperature)
      #sort the vector in ascending order
      sorted_nums <- sort(temps)
      # get the length of the vector
      n <- length(sorted_nums)

      # select the three middle values

      middle_vals <- sorted_nums[(n/2-1):(n/2 + 1)]
      df<-df[df$temperature %in% middle_vals,]
    }
  }else{
    if(any(df$Run=="OnePot")){
      return(df)
    }else{
      temps<-as.character(temps)
      df$temperature<-as.character(df$temperature)
      df<-df[df$temperature %in% temps,]
      df$temperature<-as.integer(df$temperature)
    }
  }
  df$temperature<-ifelse(stringr::str_detect(df$treatment,
                                             "[[:digit:]]+.[:digit:]"),
                                             stringr::str_extract(df$treatment,"[[:digit:]]+.[:digit:]"),
                                             stringr::str_extract(df$treatment,"[[:digit:]]+"))
  return(df)
}

