#Compute p-values from ATE using SCAM with RE
#' @importFrom stringr str_extract

compute_pvalues_ATE_RE_F = function(proteins,workers=6){

  if(any(names(proteins)=="Protein")){
  }
  if(any(names(proteins)=="Condition")&!any(names(proteins)=="treatment")){
    proteins$treatment=as.factor(proteins$Condition)
  }
  if(any(names(proteins)=="Abundance")&!any(names(proteins)=="I")){
    proteins$I=proteins$Abundance
  }
  if(any(names(proteins)=="Subject")&!any(names(proteins)=="sample_id")){
    proteins$sample_id=proteins$Subject
  }
  if(!any(names(proteins)=="sample_name")&any(names(proteins)=="Run")){
    proteins$sample_name=stringr::str_extract(unique(proteins$Run)[length(unique(proteins$Run))],"[[:lower:]]+")

  }else{
    proteins$sample_name<-"Simulation"
  }

  original_result = fit_scam_marginal_ATE_RE_F(proteins) |> dplyr::bind_rows()
  original_result$F_adjBH=stats::p.adjust(original_result$F_pvalue,method="BH")
  original_result$ATE_padjBH=stats::p.adjust(original_result$p.value,method="BH")

  #original_result$dTm<-original_result$Tm_treatment-original_result$Tm_vehicle
  return(original_result)
}
