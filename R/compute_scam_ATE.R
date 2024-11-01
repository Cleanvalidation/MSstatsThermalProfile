compute_scam_ATE = function(data){
  model = R.utils::withTimeout(scam::scam(Abundance ~ s(temperature,by=Condition,bs="mpd") + Condition + s(Subject,bs="re"),data = data,optimizer = "efs"),
                               onTimeout = "error", timeout = 20)
  ATE = marginaleffects::tidy(marginaleffects::comparisons(model,var="Condition"))
  return(ATE)
}

poss_compute_scam_ATE = purrr::possibly(.f=compute_scam_ATE,otherwise=dplyr::tibble(type="response",
                                                                                    term="Condition",contrast="treated - vehicle",
                                                                                    estimate=NA,std.error=NA,statistic=NA,p.value=NA,
                                                                                    conf.low=NA,conf.high=NA))

compute_scam_linear_ATE = function(data){
  model = R.utils::withTimeout(scam::scam(Abundance ~ s(temperature,bs="mpd") + Condition + s(Subject,bs="re"),data = data,optimizer = "efs"),
                               onTimeout = "error", timeout = 20)
  model_summary = summary(model)$p.table[2,]
  ATE = dplyr::tibble(Estimate = model_summary["Estimate"],SE=model_summary["Std. Error"],p_value = model_summary["Pr(>|t|)"])
  return(ATE)
}

poss_compute_scam_linear_ATE = purrr::possibly(.f=compute_scam_linear_ATE,otherwise=dplyr::tibble(Estimate=NA,SE=NA,p_value=NA))

compute_scam_ATE_noRE = function(data){
  model = R.utils::withTimeout(scam::scam(Abundance ~ s(temperature,by=Condition,bs="mpd") + Condition,data = data,optimizer = "efs"),
                               onTimeout = "error", timeout = 20)
  ATE = marginaleffects::tidy(marginaleffects::comparisons(model,var="Condition"))
  return(ATE)
}
