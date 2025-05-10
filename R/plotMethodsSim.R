plotMethodsSim<-function(df,labels=NA,processing="",temps=c("53.8","57.1","60.4")){
  if (!requireNamespace("MSstatsTMT", quietly = TRUE)) {
    stop("The MSstatsTMT package is required but not installed.")
  }
  if (!requireNamespace("MSstats", quietly = TRUE)) {
    stop("The MSstats package is required but not installed.")
  }
  if (!requireNamespace("MSstatsConvert", quietly = TRUE)) {
    stop("The MSstatsConvert package is required but not installed.")
  }

  if(!length(unique(df$Condition))==2){
    df$Condition<-stringr::str_extract(stringr::str_extract_all(stringr::str_to_lower(df$Condition),"_[[:lower:]]+"),"[[:lower:]]+")
    df$Condition<-ifelse(df$Condition=="character","Norm",df$Condition)
    df$Condition<-df$treatment<-ifelse(df$Condition=="vehicle","vehicle","treated")
  }else{
    df$Condition<-df$treatment<-ifelse(df$Condition=="vehicle","vehicle","treated")
    df$treatment<-df$Condition
  }
  df1<-df
  l<-labels
  num<-length(unique(df$Abundance|>na.omit()))
  df$num<-num
  if(any(is.na(labels))){
    l<-"AUTO"
  }
  df$shape<-as.factor(ifelse(df$temperature==min(df$temperature,na.rm=T),"reference","included"))
  #Start with NPARC
  NPARC_TP_Sigmoid<-NPARC::NPARCfit(x=df$temperature,
                                    y=df$Abundance,
                                    id=df$Accession,
                                    returnModels = TRUE,
                                    groupsAlt=as.character(stringr::str_extract(stringr::str_to_lower(df$Condition),"[[:lower:]]+")))


  NPARC_TP_Sigmoid$predictions<-NPARC_TP_Sigmoid$predictions |> dplyr::filter(!is.na(group))

  NPARC<-plot_NPARC_fit(df,NPARC_TP_Sigmoid,protein=df$Accession[1])+
    ggtitle(paste0(processing," proc. NPARC stat. model"))+
    ylim(0,1.5)+
    ggpubr::border(color="#FEC44F",size=2)+
    theme(text = element_text(size = 18),title = element_text(size=14))+
    xlim(min(df$temperature),70)+
    annotate(geom="text",label=paste0("n = ",num),  x=min(df$temperature,na.rm=T)+3,y=0.07,size=5)+
    annotate(geom="text",label=paste0("a. p-value: ",formatC(df$NPARC_adjpval, format = "f", digits = 3), x=min(df$temperature,na.rm=T)+8.5,y=0,size=5))
  #Then TPP
  Zebra_splinemodels<-purrr::map(list(df),function(x){
    model=lm(Abundance ~ splines::ns(temperature,df=5) + Condition,data=x)
    return(model)
  })
  TN_predictions=predict(Zebra_splinemodels[[1]],newdata=list(temperature=rep(seq(37,70,0.01),2),Condition=c(rep_len("treated",length(seq(37,70,0.01))),rep_len("vehicle",length(seq(37,70,0.01))))))

  TPP_TN_pred_df<-data.frame(temperature=rep(seq(37,70,0.01),2),Condition=as.factor(c(rep_len("treated",length(seq(37,70,0.01))),rep_len("vehicle",length(seq(37,70,0.01))))),Abundance=TN_predictions)
  TPP_TN_pred_df$shape<-as.factor(ifelse(TPP_TN_pred_df$temperature==min(TPP_TN_pred_df$temperature,na.rm=T),"reference","included"))
  TPP<-ggplot(df,mapping=aes(x=temperature,y=Abundance,color=Condition,shape=shape),size=2)+geom_point(size=2.5)+
    geom_line(data=TPP_TN_pred_df,size=1)+
    ggtitle(paste0(processing," proc. TPP stat. model"))+
    ylim(0,1.5)+
    ggpubr::border(color="#D95F0E",size=2)+
    xlim(min(df$temperature),70)+
    annotate(geom="text",label=paste0("n = ",num), x=min(df$temperature,na.rm=T)+3,y=0.07,size=5)+
    #annotate(geom="text",label=paste0("a. p-value: ",formatC(df$TPP_adjpval, format = "f", digits = 3), x=min(df$temperature,na.rm=T)+8.5,y=0,size=5))
  #Then SCAM

  df$Mixture<-df$Subject<-as.factor(df$Mixture)
  df$Condition<-as.factor(df$Condition)
  if(any(names(df)=="Subject")){
    df$Subject<-as.factor(df$Subject)
    if(length(levels(df$Subject))>length(unique(df$Subject))){
      df$Subject<-df$Mixture<-plyr::mapvalues(df$Subject,from=levels(df$Subject)[levels(df$Subject) %in% df$Subject],to=paste0("F",seq(1:length(unique(df$Subject)))))
    }
  }
  SCAMmodels<-NA
  SCAMmodels<-purrr::map(list(df),function(x){
    model=scam::scam(Abundance ~ s(temperature,by=Condition,bs="mpd",k=5)+
                       s(Mixture,bs="re") +
                       Condition,data=x,optimizer="efs")
    return(model)
  })


  SCAM_pred<-data.frame(temperature=rep(seq(37,70,0.01),2),
                        Condition=c(rep_len("treated",length(seq(37,70,0.01))),rep_len("vehicle",length(seq(37,70,0.01)))))|>
    dplyr::inner_join(data.frame(Condition=c(rep("vehicle",2),rep("treated",2)),
                                 Subject=paste0("F",seq(1,4)),
                                 Mixture=paste0("F",seq(1,4))))|>as.list()


  TP_predictions=predict(SCAMmodels[[1]],newdata=SCAM_pred)

  df2<-data.frame(SCAM_pred)|>dplyr::mutate(Abundance=TP_predictions,Subject=as.factor(Subject))|>dplyr::distinct()
  df2$shape<-as.factor(ifelse(df2$temperature==min(df2$temperature,na.rm=T),"reference","included"))
  df_SCAM<-df|>dplyr::mutate(Condition=treatment)
  df_SCAM$shape<-as.factor(ifelse(df_SCAM$temperature==min(df_SCAM$temperature,na.rm=T),"reference","included"))
  SCAM<-ggplot(df_SCAM,mapping=aes(x=temperature,y=Abundance,color=Condition,shape=shape))+
    geom_point(size=2.5)+
    geom_line(data=df2,size=1)+
    ggtitle(paste0(processing," proc. SCAM stat. model"))+
    ylim(0,1.5)+ggpubr::border(color="#07fff8",size=2)+
    xlim(min(df$temperature),70)+
    annotate(geom="text",label=paste0("n = ",num), x=min(df$temperature,na.rm=T)+3,y=0.07,size=5)+
    annotate(geom="text",label=paste0("a. p-value: ",formatC(df$SCAM_adjpval, format = "f", digits = 3), x=min(df$temperature,na.rm=T)+8.5,y=0,size=5))

  #Finally MSstatsTMT

  df<-df1
  if(!stringr::str_detect(df$Condition[1],"_")){
    df$Condition<-df$Group<-paste0(df$Channel,"_",df$Condition)
    if(!any(df$Condition=="Norm")){

      df$Condition[df$temperature==min(df$temperature)]<-"Norm"
    }

  }else{
    df$Condition<-df$Group<-paste0(df$Channel,"_",df$treatment)
  }
  # df<-df|>dplyr::group_by(Channel)|>
  #   dplyr::group_split()|>
  #   purrr::keep(function(x)sum(is.na(x$Abundance))/length(x$Abundance)<0.75)|>
  #   dplyr::bind_rows()
  Protein<-list(ProteinLevelData=df,FeatureLevelData=NA)

  comparison<-make_contrast_matrix(Protein,temps=temps)
  if(!any(names(df)=="model")){
    MSstat_TPresult = MSstatsTMT::groupComparisonTMT(data=Protein,
                                                     contrast.matrix=comparison,
                                                     moderated = FALSE,
                                                     adj.method = "BH",
                                                     remove_norm_channel = TRUE,
                                                     remove_empty_channel = FALSE,
                                                     save_fitted_models = TRUE,
                                                     use_log_file = FALSE,
                                                     append = FALSE,
                                                     verbose = TRUE,
                                                     log_file_path = NULL)$FittedModel
    name<-names(MSstat_TPresult)
  }else{
    MSstat_TPresult<-df$model
    name<-names(MSstat_TPresult)
  }
  TP_LMM<-tryCatch(MSstat_TPresult[[1]]@frame,error=function(x){
    data.frame(Abundance=predict(MSstat_TPresult[[1]]),Group=MSstat_TPresult[[1]]$model$Group)

  })

  TP_LMM$Channel<-stringr::str_remove(TP_LMM$Group,"_[[:lower:]]+")

  TP_LMM$Condition<-as.factor(stringr::str_extract(TP_LMM$Group,"[[:lower:]]+"))

  Protein<-df|>dplyr::mutate(Condition=as.factor(stringr::str_extract(Condition,"[[:lower:]]+")))


  TP_LMM<-TP_LMM|>dplyr::distinct()|>dplyr::inner_join(Protein|>dplyr::select(temperature,Channel))
  Protein$Subject<-as.factor(Protein$Subject)
  if(length(levels(Protein$Subject))>length(unique(Protein$Subject))){
    levels(Protein$Subject)<-droplevels(Protein$Subject)
  }
  if(!any(names(TP_LMM)=="Subject")&any(names(TP_LMM)=="Run")){
    TP_LMM<-TP_LMM|>dplyr::inner_join(Protein|>dplyr::select(Run,Subject))
  }else{
    TP_LMM<-TP_LMM|>
      dplyr::group_by(Abundance,Group)|>
      dplyr::mutate(Subject=unique(Protein$Subject))|>
      dplyr::ungroup()
  }

  temps<-temps[order(temps)]

  TP_LMM<-TP_LMM|>as.data.frame()|>
    dplyr::distinct()
  Protein$shape<-as.factor(ifelse(Protein$temperature==min(Protein$temperature,na.rm=T),"reference","included"))
  TP_LMM$shape<-as.factor(ifelse(TP_LMM$temperature==min(TP_LMM$temperature,na.rm=T),"reference","included"))
  Protein$Condition<-Protein$treatment
  if(max(Protein$Abundance,na.rm=T)<1.5){
    MSstat<-ggplot(Protein,mapping=aes(x=temperature,y=Abundance,color=Condition,shape=shape))+
      geom_point(size=2.5)+
      geom_step(data=TP_LMM,mapping=aes(linetype=Subject),size=1,alpha=0.5)+
      geom_rect(data=TP_LMM|>dplyr::filter(temperature %in% temps),
                aes(xmin = min(temperature),
                    xmax = max(temperature),
                    ymin = 0,
                    ymax = 1.25
                ),
                fill="#2C7FB8",
                alpha = 0.02,
                inherit.aes=FALSE)+
      ggtitle(paste0(processing," proc. MSstatsTMT stat. model"))+
      ylim(0,1.25)+
      ggpubr::border(color="#2C7FB8",size=2)+
      theme(text = element_text(size = 20),legend.position = "bottom",title = element_text(size=15))+
      xlim(min(df$temperature),70)+
      annotate(geom="text",label=paste0("n = ",num),  x=min(df$temperature,na.rm=T)+3,y=0.07,size=5)#+
      #annotate(geom="text",label=paste0("a. p-value: ",formatC(df$MSstats_adjpval, format = "f", digits = 3), x=min(df$temperature,na.rm=T)+8.5,y=0,size=5))
  }else{
    MSstat<-ggplot(Protein,mapping=aes(x=temperature,y=Abundance,color=Condition,shape=shape))+
      geom_point(size=2.5)+
      geom_step(data=TP_LMM,mapping=aes(linetype=Subject),size=1,alpha=0.5)+
      geom_rect(data=TP_LMM|>dplyr::filter(temperature %in% temps),
                aes(xmin = min(temperature),
                    xmax = max(temperature),
                    ymin = 0,
                    ymax = 1.5
                ),
                fill="#2C7FB8",
                alpha = 0.02,
                inherit.aes=FALSE)+
      ggtitle(paste0(processing," proc. MSstatsTMT stat. model"))+
      ylim(0,1.6)+
      ggpubr::border(color="#2C7FB8",size=2)+
      theme(text = element_text(size = 20),legend.position = "bottom",title = element_text(size=15))+
      xlim(min(df$temperature),70)+
      annotate(geom="text",label=paste0("n = ",num),  x=min(df$temperature,na.rm=T)+3,y=0.07,size=5)#+
      #annotate(geom="text",label=paste0("a. p-value: ",formatC(df$MSstats_adjpval, format = "f", digits = 3), x=min(df$temperature,na.rm=T)+8.5,y=0,size=5))
  }

  y<-cowplot::get_legend(MSstat)
  plotS2<-list(NPARC,TPP,SCAM,MSstat)
  if(max(TP_LMM$Abundance,na.rm=T)>1.5){
    plotS2<-lapply(plotS2,function(x) x+ylim(0,1.6))
  }else{
    plotS2<-lapply(plotS2,function(x) x+ylim(0,1.25))
  }
  plotS2<-lapply(plotS2,function(x) x+theme(text = element_text(size = 20),legend.position = "bottom",axis.title.x = element_text(size=20),axis.title.y=element_text(size=20),title = element_text(size=13)))
  P_fits<-ggpubr::ggarrange(plotlist=plotS2,ncol=4,nrow=1,font.label = list(size = 23, color = "black"),labels = l,legend.grob = y)
  return(P_fits)

}
