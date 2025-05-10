#' @importFrom stringr str_extract

make_contrast_matrix = function(data,temps=NA){

  if(!any(names(data$ProteinLevelData)=="temperature")){
    stop("No temperature column detected, please add a temperature column")
  }
  #if the condition column does not have channel and treatment, append it
  if(!any(stringr::str_detect(data$ProteinLevelData$Condition,"[:punct:]"))){
    data$ProteinLevelData$Condition<-paste0(data$ProteinLevelData$Channel,"_",data$ProteinLevelData$treatment)
  }
  if(any(names(data$ProteinLevelData)=="Condition")){
    data$ProteinLevelData$treatment[!data$ProteinLevelData$Condition=="Norm"]<-stringr::str_extract(data$ProteinLevelData$Condition[which(!data$ProteinLevelData$Condition=="Norm")],"[[:lower:]]+")
    #label the reference channel as norm
    data$ProteinLevelData$treatment[is.na(data$ProteinLevelData$treatment)]<-"Norm"
    data$ProteinLevelData<-data$ProteinLevelData[!data$ProteinLevelData$Condition=="Norm",]

  }

  data$ProteinLevelData$temperature<-as.character(data$ProteinLevelData$temperature)

  conditions_selected<-data$ProteinLevelData|>
    dplyr::select(temperature,Condition,treatment,Channel)|>
    dplyr::filter(temperature %in% temps)|>unique()
  condition_levels<-data$ProteinLevelData|>dplyr::select(temperature,Condition,treatment,Channel)
  #Channel 127C => Temp=44, Channel 130C =>Temp = 63
  null_contrasts<-condition_levels|>dplyr::mutate(ATE=rep(0,nrow(condition_levels)))
  #Check if there are two of each channel_condition
  null_contrasts<-null_contrasts|>
    dplyr::group_by(Channel)|>
    dplyr::mutate(n=dplyr::n())|>
    dplyr::filter(n>=2)|>
    dplyr::ungroup()|>dplyr::distinct()

  ATE=null_contrasts|>dplyr::mutate(ATE=ifelse(Condition %in% conditions_selected$Condition,
                                               2/nrow(conditions_selected)*ifelse(stringr::str_detect(stringr::str_to_lower(Condition)
                                                                                                      ,"treat"),-1,
                                                                                  ifelse(stringr::str_detect(stringr::str_to_lower(Condition)
                                                                                                             ,"vehicle"),1,ATE)),
                                               ATE))|>dplyr::select(ATE)
  #Treated - Vehicle
  #levels_Condition = unique(data$ProteinLevelData$Group)

  contrast_matrix = t(as.matrix(ATE))
  if(length(temps)>=10){
    colnames(contrast_matrix) = unique(null_contrasts$Condition)
  }else{
    names(contrast_matrix)= null_contrasts$Condition
  }
  rownames(contrast_matrix) = c("ATE")
  colnames(contrast_matrix)<-unique(null_contrasts$Condition)
  return(contrast_matrix)
}
