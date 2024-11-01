#Compute p-values from ATE using SCAM with RE
compute_pvalues_ATE_RE_F = function(fulldata,workers=6){
  if(any(names(fulldata)=="Protein")){
    fulldata$Accession=fulldata$Protein
  }
  if(any(names(fulldata)=="Condition")&!any(names(fulldata)=="treatment")){
    fulldata$treatment=as.factor(fulldata$Condition)
  }
  if(any(names(fulldata)=="Abundance")&!any(names(fulldata)=="I")){
    fulldata$I=fulldata$Abundance
  }
  if(any(names(fulldata)=="Subject")&!any(names(fulldata)=="sample_id")){
    fulldata$sample_id=fulldata$Subject
  }
  if(!any(names(fulldata)=="sample_name")&any(names(fulldata)=="Run")){
    fulldata$sample_name=stringr::str_extract(unique(fulldata$Run)[length(unique(fulldata$Run))],"[[:lower:]]+")

  }else{
    fulldata$sample_name<-"Simulation"
  }

  original_result = fit_scam_marginal_ATE_RE_F(fulldata) |> dplyr::bind_rows()
  original_result$F_adjBH=stats::p.adjust(original_result$F_pvalue,method="BH")
  original_result$ATE_padjBH=stats::p.adjust(original_result$p.value,method="BH")

  #original_result$dTm<-original_result$Tm_treatment-original_result$Tm_vehicle
  return(original_result)
}
