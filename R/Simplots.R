#' @importFrom stringr str_extract

norm_exp_data_sigmoid<-function(x){
  df1<-Human_O00267|>dplyr::mutate(Protein=Accession)
  df1$shape<-as.factor(ifelse(df1$temperature==min(df1$temperature,na.rm=T),"reference","included"))
  #Start with NPARC
  NPARC_TN_Sigmoid<-NPARC::NPARCfit(x=df1$temperature,
                                    y=df1$Abundance,
                                    id=df1$Accession,
                                    returnModels = TRUE,
                                    groupsAlt=as.character(stringr::str_extract(stringr::str_to_lower(df1$Condition),"[[:lower:]]+")))


  NPARC_TN_Sigmoid$predictions<-NPARC_TN_Sigmoid$predictions |> dplyr::filter(!is.na(group))
  NPARC_theo<-NPARC::NPARCtest(NPARC_TN_Sigmoid$metrics,dfType="theoretical")

  NPARC<-plot_NPARC_fit(df1,NPARC_TN_Sigmoid,protein=df$Accession[1])+
    ggtitle(paste0(processing," proc. NPARC stat. model"))+
    ylim(0,1.5)+
    ggpubr::border(color="#FEC44F",size=2)+
    theme(text = element_text(size = 18),title = element_text(size=14))+
    xlim(min(df$temperature),70)+
    annotate(geom="text",label=paste0("n = ",num),  x=min(df$temperature,na.rm=T)+3,y=0.07,size=5)+
    annotate(geom="text",label=paste0("p-value: ",formatC(NPARC_theo$pVal, format = "e", digits = 2)), x=min(df$temperature,na.rm=T)+6.3,y=0,size=5)
  #generate a model


  #Then TPP
  Zebra_splinemodels<-purrr::map(list(df1),function(x){
    model=lm(Abundance ~ splines::ns(temperature,df=5) + Condition,data=x)
    return(model)
  })
  TN_predictions=predict(Zebra_splinemodels[[1]],newdata=list(temperature=rep(seq(37,70,0.01),2),Condition=c(rep_len("treated",length(seq(37,70,0.01))),rep_len("vehicle",length(seq(37,70,0.01))))))

  TPP_TN_pred_df<-data.frame(temperature=rep(seq(37,70,0.01),2),Condition=as.factor(c(rep_len("treatment",length(seq(37,70,0.01))),rep_len("vehicle",length(seq(37,70,0.01))))),Abundance=TN_predictions)
  TPP_TN_pred_df$shape<-as.factor(ifelse(TPP_TN_pred_df$temperature==min(TPP_TN_pred_df$temperature,na.rm=T),"reference","included"))

  df<-TPPnorm_Human_Proteins
  df$Condition<-ifelse(stringr::str_detect(stringr::str_to_lower(df$Experiment),"vehicle"),
                                                                     "vehicle","treatment")
  df$shape<-as.factor(ifelse(df$temperature==min(df$temperature,na.rm=T),"reference","included"))
  df$Condition<-paste0(df$Channel,"_",df$Condition)
  TPP_trueneg_results = TPP_NPARC_calc(df,method="NPARC",DF=5,temps=set_temps(10,unique(df$temperature)),
                                       CARRIER=FALSE,returnModels=TRUE,filters=FALSE,NORM = FALSE)
  TPP_trueneg_results<-TPP_trueneg_results|>dplyr::filter(Protein %in% df1$Protein[1])
  df<-df|>dplyr::filter(Protein %in% df1$Protein[1])
  df$Condition<-as.factor(ifelse(stringr::str_detect(stringr::str_to_lower(df$Condition),"vehicle"),
                                 "vehicle","treatment"))
  TPP<-ggplot(df,mapping=aes(x=temperature,y=Abundance,color=Condition,shape=shape),size=2)+geom_point(size=2.5)+
    geom_line(data=TPP_TN_pred_df,size=1)+
    ggtitle(paste0(processing," proc. TPP stat. model"))+
    ylim(0,1.5)+
    ggpubr::border(color="#D95F0E",size=2)+
    xlim(min(df$temperature),70)+
    annotate(geom="text",label=paste0("n = ",num), x=min(df$temperature,na.rm=T)+3,y=0.07,size=5)+
    annotate(geom="text",label=paste0("p-value: ",formatC(TPP_trueneg_results$pNPARC[1], format = "e", digits = 2)), x=min(df$temperature,na.rm=T)+6.3,y=0,size=5)
}
