#'
#' fit_scam_marginal_ATE_RE_F
#'
#'@description Fits Shape constrained additive models to calculate average treatment effect
#'@importFrom stringr str_extract
#'@importFrom marginaleffects comparison


fit_scam_marginal_ATE_RE_F= function(accession_data){

  if(any(names(accession_data)=="Mixture")&any(names(accession_data)=="sample_name")){
    if(!any(names(accession_data)=="treatment")&any(names(accession_data)=="Condition")){
      accession_data$Condition<-accession_data$treatment<-factor(stringr::str_extract(stringr::str_to_lower(accession_data$treatment),"[[:lower:]]+"),levels=c("vehicle","treated"))
    }else if(any(stringr::str_detect(accession_data$treatment[1],"_"))){
      accession_data$Condition<-accession_data$treatment<-factor(stringr::str_extract(stringr::str_to_lower(accession_data$treatment),"[[:lower:]]+"),levels=c("vehicle","treated"))
    }
    if(!any(names(accession_data)=="Protein")){
      accession_data$Accession<-accession_data$Protein
    }
    accession_data<-accession_data |> dplyr::select(Accession,I,temperature,sample_name,treatment,sample_id,Condition)
    data_gdf_accession = accession_data |> dplyr::group_split(Accession,sample_name)
  }else{
    data_gdf_accession = accession_data |> dplyr::group_split(Accession)

  }
  model<-furrr::future_map(data_gdf_accession,function(x) tryCatch(poss_fit_scam_ANOVA_RE(x),error=function(cond){error_df = tibble(type=NA,
                                                                                                                                    term=NA,
                                                                                                                                    estimate=NA,
                                                                                                                                    std.error=NA,
                                                                                                                                    statistic=NA,
                                                                                                                                    conf.low=NA,
                                                                                                                                    conf.high=NA,
                                                                                                                                    adj_r2=NA,
                                                                                                                                    Accession=accession_data$Accession[1],
                                                                                                                                    slope.pval=NA,
                                                                                                                                    Tm_treated=NA,
                                                                                                                                    Tm_vehicle=NA
  )
  return(error_df)}
  ))
  model<-model %>% purrr::keep(function(x) length(x)>1)
  model1<-purrr::map(model,function(x) x[1])
  model2<-purrr::map(model,function(x) x[2])#this is the actual model


  model2<-purrr::map(model2,function(x) x[[1]])
  model2<-model2#%>% purrr::keep(function(x) any(class(x)=="scam"))
  ATE<-model2|>purrr::keep(function(x) any(class(x)=="scam"))

  check_ATE<-furrr::future_map(model2,function(x)tryCatch(broom::tidy(marginaleffects::comparisons(x,variables="treatment")),
                                               error=function(e){return(FALSE)}))
  model1<-purrr::map(model1,function(x) x[[1]])
  check<-purrr::map2(model1,check_ATE,function(x,y){tryCatch(
    base::cbind(x,y),error=function(e){return(data.frame())})})



  ATE=check

  return(ATE)

}
