make_profiles<-function(dataList,n_protein,s_Var,r_Var){
  stopifnot(length(r_Var)==1|!length(s_Var==5))
  dataList$uniqueID<-dataList$Protein

  # Generate a vector of subject IDs
  subject <- rep(paste0("F",seq(1:4)), each = 10)
  TechRepMixture<-data.frame(Subject=c(rep("Vehicle",2),rep("Treatment",2)))|>
    dplyr::group_by(Subject)|>
    dplyr::mutate(TechRepMixture=sample(c(1,2),size=2,replace=FALSE))|>
    dplyr::ungroup()|>
    dplyr::mutate(Condition=Subject,
                  Subject=paste0(Subject,"_",TechRepMixture))

  #use the first protein fitted values on the list as a template
  fit<-dataList$fit[[1]]
  #transform fitted values to vectors
  fitted_values<-data.frame(Abundance=t(predict(fit))|>
                              as.vector())
  fitted_values$temperature<-NA
  if(length(fitted_values$Abundance)==20){
    # Generate a vector of temps
    temps <- unique(dataList$temperature)[1:10]
    fitted_values$temperature<-rep(as.character(temps),2)
  }else{
    # Generate a vector of temps
    temps <- unique(dataList$temperature)[2:10]
    fitted_values$temperature<-rep(as.character(temps[1:(length(fitted_values$Abundance)/2)]),2)
    if(!length(fitted_values$temperature)==20){
      temps<-data.frame(temperature=rep(as.character(unique(dataList$temperature)[1:10]),2))
      fitted_values$temperature<-as.character(fitted_values$temperature)
      fitted_values<-fitted_values|>dplyr::right_join(temps)
      fitted_values$Abundance<-ifelse(fitted_values$temperature==min(unique(dataList$temperature),na.rm=T),1,fitted_values$Abundance)
    }
  }

  # Create a data frame with the temps, subject IDs, group labels, and fitted values
  df <- data.frame(temperature =rep(temps$temperature,4),
                   Condition=c(rep("Vehicle",20),rep("Treatment",each=20)))|>
    dplyr::group_by(temperature,Condition)|>
    dplyr::inner_join(TechRepMixture)|>
    dplyr::ungroup()

  df$temperature<-as.character(df$temperature)

  df<-rep(list(data.frame(df)),length(unique(dataList$Protein)))
  names(df)<-unique(dataList$Protein)

#generate fitted values for the selected protein
  fitted_values<-rep(list(data.frame(fitted_values)),length(unique(dataList$Protein)))
  names(fitted_values)<-unique(dataList$Protein)

  df1<-purrr::map2(df,fitted_values,function(x,y){
    z<-x|>dplyr::left_join(y,by="temperature")|>
      dplyr::distinct()|>
      dplyr::mutate(Abundance=ifelse(is.na(Abundance),1,Abundance))|>
      dplyr::distinct()|>
      dplyr::group_by(Subject,Condition)|>
      dplyr::filter(!duplicated(temperature))
    return(z)
  })
  data<-list(unique(dataList$Protein))
  df<-purrr::map2(df1,data,function(x,y)x|>dplyr::mutate(uniqueID=y))

  df<-dplyr::bind_rows(df)

  #partition data by temperature to ensure 2 bioreps for vehicle and 2 bioreps for treated
  df<-df|>
    dplyr::mutate(Condition=as.factor(Condition),
                  temperature=as.numeric(temperature))|>
    dplyr::distinct()

  #how many times will the proteins be simulated
  fitted_list<-rep(list(data.frame(df)),n_protein)
  fitted_listnames<-data.frame(uniqueID=paste0("Sim_",seq(1:n_protein)))|>dplyr::group_by(uniqueID)|>
    dplyr::group_split()

  temps<-dataList|>dplyr::select(temperature,Channel)|>dplyr::distinct()

  for(z in seq(s_Var)){
    Sims<-purrr::map2(fitted_list,fitted_listnames,function(x,y)add_variation_shift(x,icc = s_Var[z]/(r_Var+s_Var[z]),re_var=r_Var,t_range = temps$temperature)|>
                        dplyr::mutate(Protein=y$uniqueID,uniqueID=y$uniqueID,bio_var = s_Var[z],re_var=r_Var,icc=s_Var[z]/(r_Var+s_Var[z])))
  }
  Sims<-dplyr::bind_rows(Sims)|>dplyr::inner_join(temps)
  Sims$Condition<-as.character(Sims$Condition)
  return(Sims)
}
