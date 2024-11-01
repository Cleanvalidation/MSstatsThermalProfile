#This function fits a scam function with random effect on the subject
fit_scam_ANOVA_RE = function(accession_data){
  accession_data<-accession_data %>% dplyr::mutate(sample_id=as.factor(sample_id))
  if(any(!class(accession_data$treatment)=="factor")){
  accession_data$treatment<-as.factor(accession_data$treatment)
  accession_data$sample_id<-as.factor(accession_data$sample_id)
  }
  mod1 = R.utils::withTimeout(scam::scam(I ~ s(temperature,by=treatment,bs="mpd",k=5)+
                                           s(sample_id,bs="re") +
                                           treatment,data=accession_data,optimizer="efs"),timeout=120)
  mod = R.utils::withTimeout(scam::scam(I ~ s(temperature,bs="mpd",k=5)+
                                          s(sample_id,bs="re"),data=accession_data,optimizer="efs"),timeout=120)
  if(any(class(mod1)=="scam") &any(class(mod)=="scam")){
    anova_F<-R.utils::withTimeout(anova(mod,mod1,test="F"),timeout=120)
    anova_F<-anova_F[2,] %>% as.data.frame()
    anova_F<-data.frame(Accession=accession_data$Accession[1],r.sq=summary(mod1)$r.sq,dRSS=anova_F$Deviance,F_test=anova_F$F,F_pvalue=anova_F$`Pr(>F)`)
    anova_F<-list(anova_F,mod1)

  }else{
    anova_F<-mod1
  }

  return(anova_F)
}
poss_fit_scam_ANOVA_RE = purrr::possibly(.f=fit_scam_ANOVA_RE,otherwise="Error")
