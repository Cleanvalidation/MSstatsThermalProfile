
make_contrast_matrix_all = function(data,temps=NA){
  if(!any(names(data$ProteinLevelData)=="temperature")){
    stop("No temperature column detected, please add a temperature column")
  }
  #data<-list(ProteinLevelData=choose_temps(data$ProteinLevelData))

  if(any(stringr::str_detect(data$ProteinLevelData$Condition,"Norm"))){
    #remove norm condition
    data$ProteinLevelData<-data$ProteinLevelData[!data$ProteinLevelData$Condition=="Norm",]
  }
  # if("131C" %in% unique(data$ProteinLevelData$Channel)){
  #
  #   stop("Remove Reference Channel 131C")
  # }
  contrasts<-choose_temps(data$ProteinLevelData,temps=unique(data$ProteinLevelData$temperature))
  if(!is.na(temps)){
    condition_levels<-contrasts|>dplyr::select(temperature,Condition,Mixture)|>
      unique()
  }else{
    condition_levels<-contrasts|>dplyr::select(temperature,Condition)|>unique()
  }
  #Channel 127C => Temp=44, Channel 130C =>Temp = 63
  null_contrasts<-data.frame(Condition=unique(condition_levels$Condition),
                             ATE=rep(0,nrow(data$ProteinLevelData|>dplyr::select(temperature,Condition)|>unique())))|>dplyr::distinct()

  ATE=null_contrasts|>dplyr::mutate(ATE=ifelse(Condition %in% unique(condition_levels$Condition),
                                               2/nrow(unique(condition_levels))*ifelse(stringr::str_detect(stringr::str_to_lower(Condition)
                                                                                                   ,"vehicle"),-1,1),
                                               ATE))|>dplyr::select(ATE)
  #Treated - Vehicle
  #levels_Condition = unique(data$ProteinLevelData$Group)

  contrast_matrix = t(as.matrix(ATE))


  col_names<-unique(null_contrasts$Condition)

  colnames(contrast_matrix) = col_names

  rownames(contrast_matrix) = c("ATE")
  return(contrast_matrix)
}
