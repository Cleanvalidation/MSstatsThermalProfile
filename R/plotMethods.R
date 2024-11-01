plotMethods<-function(target_protein,labels=NA,processing="",temps=c("53.8","57.1","60.4"),fit="Spline"){
  #since the MSstats protein summarization output does not have the temperature data, we need to perform an inner join to add this
  target_protein$treatment<-stringr::str_extract(stringr::str_to_lower(target_protein$treatment),"[[:lower:]]+")
  target_protein$treatment<-ifelse(stringr::str_detect(target_protein$treatment,"ehicle"),"vehicle","treated")


  target_protein<-dplyr::bind_rows(target_protein)|>
    dplyr::mutate(Condition=as.character(Condition),
                  Accession=as.character(Protein),
                  Protein=as.character(Protein),
                  Abundance=as.numeric(Abundance))
  target_protein$temperature<-as.numeric(target_protein$temperature)

  target_protein1<-target_protein

  target_protein1<-target_protein1|>
    dplyr::mutate(shape=ifelse(temperature==min(temperature,na.rm=TRUE),"reference","included"))

  if(processing=="MSstatsTMT"){#If this is MSstatsTMT processed, unlog and ratio to reference channel is needed for benchmark
    target_protein1$Condition<-target_protein1$treatment
    target_protein1$Experiment<-paste0(target_protein1$treatment,"_",target_protein1$TechRepMixture)

    target_protein1<-target_protein1|>
      dplyr::mutate(Abundance=2^Abundance)|>
      dplyr::bind_rows()

    target_protein1$Abundance<-suppressWarnings(target_protein1$Abundance/target_protein1$Abundance[target_protein1$shape=="reference"])

  }else if(any(names(target_protein1)=="treatment")&any(names(target_protein1)=="Experiment")){#if this is TPP processed, no ratio for the reference channel is needed
    target_protein1$Condition<-stringr::str_extract(stringr::str_to_lower(as.character(target_protein1$Experiment)),"[[:lower:]]+")
    target_protein1$treatment<-target_protein1$Condition

  }else if(any(names(target_protein1)=="Experiment")&any(stringr::str_detect(target_protein1$Condition[1],"_"))){#if this is TPP processed without a treatment column
    target_protein1$Condition<-stringr::str_extract(stringr::str_to_lower(as.character(target_protein1$Experiment)),"[[:lower:]]+")
    target_protein1$treatment<-target_protein1$Condition

  }



  l<-labels
  num<-length(unique(target_protein1$Abundance|>na.omit()))
  target_protein1$num<-num
  if(any(is.na(labels))){
    l<-"AUTO"
  }
  target_protein1$shape<-as.factor(ifelse(target_protein1$temperature==min(target_protein1$temperature,na.rm=T),"reference","included"))
  if(!any(names(target_protein1)=="Accession")&any(names(target_protein1)=="uniqueID")){
    target_protein1$Accession<-target_protein1$uniqueID
  }

  #Start with NPARC
  NPARC_Sigmoid<-NPARC::NPARCfit(x=target_protein1$temperature,
                                 y=target_protein1$Abundance,
                                 id=target_protein1$Accession,
                                 returnModels = TRUE,
                                 groupsAlt=target_protein1$treatment)


  NPARC_Sigmoid$predictions<-NPARC_Sigmoid$predictions |> dplyr::filter(!is.na(group))
  target_protein1$NPARC_adjpval<-as.numeric(target_protein1$NPARC_adjpval)

  NPARC<-plot_NPARC_fit(target_protein1,NPARC_Sigmoid,protein=target_protein1$Accession[1])+
    ggtitle(paste0(processing," proc. NPARC stat. model"))+
    ylim(0,1.5)+
    ggpubr::border(color="#FEC44F",size=2)+
    annotate(geom="text",label=paste0("a. p-value: ",sprintf("%.3f",target_protein1$NPARC_adjpval[1])), x=min(target_protein1$temperature,na.rm=T)+5.3,y=0,size=5)

  #Then TPP
  if(length(unique(target_protein1$Condition)==1)){
    target_protein1$Condition<-as.character(target_protein1$Condition)
    target_protein1$Condition<-target_protein1$treatment
  }

  if(fit=="Spline"){
    TPP_splinemodels<-purrr::map(list(target_protein1),function(x){
      model=lm(Abundance ~ splines::ns(temperature,df=5) + Condition,data=x)
      return(model)
    })
    target_protein_list<-list(target_protein1)

    new_data<-target_protein1|>dplyr::select(Condition,Experiment)|>dplyr::distinct()
    # Generate temperatures for each row
    temperatures <- seq(37, 67, length.out = 50)

    new_data<-new_data|>dplyr::group_by(Experiment)|>dplyr::reframe(temperature=temperatures,Condition=Condition)

    TPP_pred=predict(TPP_splinemodels[[1]],newdata=new_data)

    TPP_predicted_df<-new_data|>dplyr::mutate(Abundance=TPP_pred)


    target_protein1$shape<-as.factor(ifelse(target_protein1$temperature==min(target_protein1$temperature,na.rm=T),"reference","included"))
    TPP_predicted_df$shape<-as.factor(ifelse(TPP_predicted_df$temperature==min(TPP_predicted_df$temperature,na.rm=T),"reference","included"))
    TPP_predicted_df$Condition<-as.factor(TPP_predicted_df$Condition)
    target_protein1$Condition<-as.factor(target_protein1$Condition)
    TPP_predicted_df<-TPP_predicted_df|>
      dplyr::arrange(Condition, temperature)
    #Set p-value as numeric
    target_protein1$TPP_adjpval<-as.numeric(target_protein1$TPP_adjpval)

    TPP<-ggplot(target_protein1,mapping=aes(x=temperature,y=Abundance,shape=shape,color=Condition,group=Experiment))+
      geom_point()+
      geom_line(data=TPP_predicted_df
                ,linewidth=1,mapping=aes(x=temperature,y=Abundance,color=Condition))+
      ggtitle(paste0(processing," proc. TPP stat. model"))+
      ylim(0,1.5)+
      ggpubr::border(color="#D95F0E",size=2)+
      xlim(min(target_protein1$temperature),70)+
    annotate(geom="text",label=paste0("a. p-value: ",sprintf("%.3f",target_protein1$TPP_adjpval[1])), x=min(target_protein1$temperature,na.rm=T)+5.3,y=0,size=5)
  }else if(fit=="Sigmoid"){
    #partition proteins by replicate
    target_protein_list<-target_protein1|>
      dplyr::group_by(Experiment)|>
      dplyr::mutate(Condition=as.factor(treatment))|>
      dplyr::group_split()

    TPPfit<-lapply(target_protein_list,function(x) fitSigmoidTR(xVec=x$temperature,x$Abundance,startPars = c("Pl"=0, "a"=550, "b"=10),maxAttempts=50))
    x_pred1<-seq(37,70,length.out=100)
    TPPpred<-lapply(TPPfit,function(x) predict(x, newdata=list(x=x_pred1)))
    TPP_predicted_df<-purrr::map2(target_protein_list,TPPpred,function(x,y)data.frame(Abundance=y,
                                                                                      temperature=seq(37,70,length.out=100))|>
                                    dplyr::mutate(Experiment=x$Experiment[1]))
    TPP_pred_target_protein1<-dplyr::bind_rows(TPP_predicted_df)
    target_protein1$shape<-as.factor(ifelse(target_protein1$temperature==min(target_protein1$temperature,na.rm=T),"reference","included"))
    TPP_pred_target_protein1$shape<-as.factor(ifelse(TPP_pred_target_protein1$temperature==min(TPP_pred_target_protein1$temperature,na.rm=T),"reference","included"))

    TPP_pred_target_protein1$Experiment<-as.factor(TPP_pred_target_protein1$Experiment)
    TPP_pred_target_protein1$Condition<-TPP_pred_target_protein1$Experiment
    target_protein1$Condition<-as.factor(target_protein1$Condition)


    TPP<-ggplot2::ggplot(target_protein1,mapping=aes(x=temperature,y=Abundance,shape=shape,group=Experiment,color=Condition),size=2)+geom_point(size=2.5)+
      geom_line(data=TPP_pred_target_protein1,linewidth=1)+
      ggtitle(paste0(processing," proc. TPP stat. model"))+
      ylim(0,1.5)+
      ggpubr::border(color="#D95F0E",size=2)

  }

  #Then SCAM

  target_protein1$Mixture<-target_protein1$Subject<-as.factor(target_protein1$Mixture)
  target_protein1$Condition<-as.factor(target_protein1$Condition)
  if(any(names(target_protein1)=="Subject")){
    target_protein1$Subject<-as.factor(target_protein1$Subject)
    if(length(levels(target_protein1$Subject))>length(unique(target_protein1$Subject))){
      target_protein1$Subject<-target_protein1$Mixture<-plyr::mapvalues(target_protein1$Subject,from=levels(target_protein1$Subject)[levels(target_protein1$Subject) %in% target_protein1$Subject],to=paste0("F",seq(1:length(unique(target_protein1$Subject)))))
    }
  }

  SCAMmodels<-NA
  SCAMmodels<-purrr::map(list(target_protein1),function(x){
    model=scam::scam(Abundance ~ s(temperature,by=Condition,bs="mpd",k=5)+
                       s(Mixture,bs="re") +
                       Condition,data=x,optimizer="efs")
    return(model)
  })


  SCAM_pred<-data.frame(temperature=rep(seq(37,70,0.01),2),
                        Condition=c(rep_len("treated",length(seq(37,70,0.01))),rep_len("vehicle",length(seq(37,70,0.01)))))|>
    dplyr::inner_join(data.frame(Condition=c(rep("vehicle",2),rep("treated",2)),
                                 Subject=paste0("F",seq(1,4)),
                                 Mixture=paste0("F",seq(1,4)))|>
                        dplyr::mutate(treatment=Condition),relationship="many-to-many")|>as.list()


  SCAM_predictions=predict(SCAMmodels[[1]],newdata=SCAM_pred)

  target_protein2<-data.frame(SCAM_pred)|>dplyr::mutate(Abundance=SCAM_predictions,Subject=as.factor(Subject))|>dplyr::distinct()
  target_protein2$shape<-as.factor(ifelse(target_protein2$temperature==min(target_protein2$temperature,na.rm=T),"reference","included"))
  target_protein_SCAM<-target_protein1
  target_protein_SCAM$shape<-as.factor(ifelse(target_protein_SCAM$temperature==min(target_protein_SCAM$temperature,na.rm=T),"reference","included"))
  #Set p-value as numeric
  target_protein_SCAM$SCAM_adjpval<-as.numeric(target_protein_SCAM$SCAM_adjpval)
  SCAM<-ggplot(target_protein_SCAM,mapping=aes(x=temperature,y=Abundance,color=Condition,shape=shape))+
    geom_point(size=2.5)+
    geom_line(data=target_protein2,linewidth=1)+
    ggtitle(paste0(processing," proc. SCAM stat. model"))+
    ylim(0,1.5)+ggpubr::border(color="#07fff8",size=2)+
    annotate(geom="text",label=paste0("a. p-value: ",sprintf("%.3f",target_protein_SCAM$SCAM_adjpval)), x=min(target_protein$temperature,na.rm=T)+5.3,y=0,size=5)

  #Finally MSstatsTMT
  if(any(names(target_protein)=="Experiment")){
    target_protein$Condition<-paste0(target_protein$Channel,"_",stringr::str_extract(target_protein$Experiment,"[[:lower:]]+"))

    target_protein$BioReplicate<-target_protein$Condition
    target_protein$treatment<-stringr::str_extract(target_protein$Experiment,"[[:lower:]]+")

  }else{
    target_protein$Condition<-paste0(target_protein$Channel,"_",stringr::str_extract(target_protein$treatment,"[[:lower:]]+"))
  }
  #Set reference channel
  target_protein$Condition<-ifelse(target_protein$temperature==min(target_protein$temperature),"Norm",target_protein$Condition)
  #filter out reference channel from data
  Protein<-list(ProteinLevelData=target_protein[!target_protein$Condition=="Norm",],FeatureLevelData=NA)

  comparison<-make_contrast_matrix(Protein,temps=temps)
  if(!any(names(target_protein)=="model")){
    MSstat_result = groupComparisonThermalProfiling(
      Protein,
      contrast.matrix = comparison,
      missing_timepoint = "replace",
      replacement = list("131_vehicle" = c("130C_vehicle","130N_vehicle")),
      moderated = FALSE,
      adj.method = "BH",
      remove_norm_channel = TRUE,
      remove_empty_channel = FALSE,
      save_fitted_models = TRUE,
      use_log_file = FALSE,
      append = FALSE,
      verbose = TRUE,
      log_file_path = NULL
    )$FittedModel
    name<-names(MSstat_result)
  }else{
    MSstat_result<-target_protein$model
    name<-names(MSstat_result)
  }
  LMM_MSstat<-tryCatch(MSstat_result[[1]]@frame,error=function(x){
    data.frame(Abundance=predict(MSstat_result[[1]]),Group=MSstat_result[[1]]$model$Group)

  })

  LMM_MSstat$Channel<-stringr::str_remove(LMM_MSstat$Group,"_[[:lower:]]+")

  LMM_MSstat$Condition<-stringr::str_remove(as.factor(stringr::str_extract(LMM_MSstat$Group,"_[[:lower:]]+")),"_")
  if(any(names(target_protein)=="treatment")){
    Protein<-target_protein|>dplyr::mutate(Condition=as.factor(treatment),
                                           Subject=Mixture)
  }else if(any(names(target_protein)=="Experiment")){
    Protein<-target_protein|>dplyr::mutate(Condition=stringr::str_extract(as.character(Experiment),"[[:lower:]]+"),
                                           Subject=Mixture)
    Protein$treatment<-Protein$Condition
  }

  if(length(levels(Protein$Subject))>length(unique(Protein$Subject))){
    Protein$Subject<-as.factor(Protein$Subject)
    levels(Protein$Subject)<-droplevels(Protein$Subject)
  }
  temps<-temps[order(temps)]

  LMM_MSstat<-LMM_MSstat|>as.data.frame()|>
    dplyr::distinct()|>dplyr::inner_join(Protein |>dplyr::select(temperature,Channel),relationship="many-to-many")
  Protein$shape<-as.factor(ifelse(Protein$temperature==min(Protein$temperature,na.rm=T),"reference","included"))
  LMM_MSstat$shape<-"included"
  if(length(unique(Protein$Condition))>2|length(unique(Protein$Condition))==1){

    Protein$Condition<-ifelse(stringr::str_detect(
      stringr::str_extract(
        stringr::str_remove_all(
          stringr::str_to_lower(Protein$Run),"[[:digit:]]+&[[:punct:]]+"),"[[:lower:]]+")
      ,"dmso|control"),"vehicle","treated")
  }
  if(!any(names(LMM_MSstat)=="Subject")){
    LMM_MSstat$Subject<-as.factor(LMM_MSstat$Run)
  }

  #set p-value as numeric
  target_protein$MSstats_adjpval<-as.numeric(target_protein$MSstats_adjpval)

  if(max(Protein$Abundance,na.rm=T)<1.5){
    MSstat<-ggplot(Protein,mapping=aes(x=temperature,y=Abundance,color=Condition,shape=shape))+
      geom_point(aes(size=shape,color=Condition))+
      geom_step(data=LMM_MSstat,mapping=aes(linetype=Subject),size=1,alpha=0.5)+
      annotate("rect", xmin=min(as.numeric(temps)), xmax = max(as.numeric(temps)), ymin=0,ymax = 1.25, alpha=0.2, fill="#2C7FB8")+
      ggtitle(paste0(processing," proc. MSstatsTMT stat. model"))+
      ylim(0,1.25)+
      ggpubr::border(color="#2C7FB8",size=2)+
      theme(legend.position = "bottom")+
      annotate(geom="text",label=paste0("a. p-value: ",sprintf("%.3f",target_protein$MSstats_adjpval)), x=min(target_protein$temperature,na.rm=T)+5.3,y=0,size=5)
  }else if (max(Protein$Abundance,na.rm=T)>10&max(target_protein$Abundance,na.rm=T)<1.26){
    MSstat<-ggplot(Protein,mapping=aes(x=temperature,y=Abundance,color=Condition,shape=shape))+
      ylab(expression(log[2]~Abundance))+
      geom_point(aes(size=shape,color=Condition))+
      geom_step(data=LMM_MSstat,mapping=aes(linetype=Subject),linewidth=1,alpha=0.5)+
      annotate("rect", xmin=min(as.numeric(temps),na.rm=T), xmax = max(as.numeric(temps),na.rm=T), ymin=0,ymax=max(target_protein$Abundance)+0.5, alpha=0.2, fill="#2C7FB8")+
      ggtitle(paste0(processing," proc. MSstatsTMT stat. model"))+
      ggpubr::border(color="#2C7FB8",size=2)+
      theme(legend.position = "bottom")+

      ylim(min(target_protein$Abundance),max(target_protein$Abundance,na.rm=T)+1)+
      annotate(geom="text",label=paste0("a. p-value: ",sprintf("%.3f",target_protein$MSstats_adjpval)),
               x=min(target_protein$temperature,na.rm=T)+5.3,
               y=min(target_protein$Abundance,na.rm=T),size=5)
  }else if(max(target_protein$Abundance,na.rm=T)>1.4&max(Protein$Abundance,na.rm=T)>10){

    MSstat<-ggplot(Protein,mapping=aes(x=temperature,y=Abundance,color=Condition,shape=shape))+
      ylab(expression(log[2]~Abundance))+
      geom_point(aes(size=shape,color=Condition))+
      geom_step(data=LMM_MSstat,mapping=aes(linetype=Subject),linewidth=1,alpha=0.5)+
      annotate("rect",
               xmin=min(as.numeric(temps)),
               xmax = max(as.numeric(temps)),
               ymin=min(target_protein$Abundance,na.rm=T),
               ymax=max(target_protein$Abundance,na.rm=T)+0.5, alpha=0.2, fill="#2C7FB8")+
      ggtitle(paste0(processing," proc. MSstatsTMT stat. model"))+
      ylim(0,max(target_protein$Abundance,na.rm=T))+
      ggpubr::border(color="#2C7FB8",size=2)+
      theme(legend.position = "bottom")+
      ylim(min(target_protein$Abundance),max(target_protein$Abundance)+1)+
      annotate(geom="text",label=paste0("a. p-value: ",sprintf("%.3f",target_protein$MSstats_adjpval)),
               x=min(target_protein$temperature,na.rm=T)+8.3,
               y=min(target_protein$Abundance,na.rm=T),size=5)


  }else{
    MSstat<-ggplot(Protein,mapping=aes(x=temperature,y=Abundance,color=Condition,shape=shape))+
      geom_point(aes(size=shape,color=Condition))+
      geom_step(data=LMM_MSstat,mapping=aes(linetype=Subject),size=1,alpha=0.5)+
      annotate("rect", xmin=min(as.numeric(temps)), xmax = max(as.numeric(temps)), ymin=min(Protein$Abundance),ymax=max(Protein$Abundance,na.rm=T), alpha=0.2, fill="#2C7FB8")+
      ggtitle(paste0(processing," proc. MSstatsTMT stat. model"))+
      ylim(min(target_protein$Abundance,na.rm=T),max(target_protein$Abundance,na.rm=T))+
      ylab(expression(log[2]~Abundance))+
      ggpubr::border(color="#2C7FB8",size=2)+
      theme(legend.position = "bottom")+
      annotate(geom="text",label=paste0("a. p-value: ",sprintf("%.3f",target_protein$MSstats_adjpval)),
               x=min(target_protein$temperature,na.rm=T)+5.3,
               y=min(target_protein$Abundance,na.rm=T),size=5)
  }

  y<-cowplot::get_legend(MSstat)
  plotS2<-list(TPP,NPARC,SCAM,MSstat)
  if(!any(stringr::str_detect(unique(Protein$temperature),"67"))){
    temps<-as.numeric(c(unique(Protein$temperature),"67"))
  }else{
    temps<-unique(Protein$temperature)
  }
  if(max(LMM_MSstat$Abundance,na.rm=T)>1.5&max(LMM_MSstat$Abundance,na.rm=10)<10){
    plotS2[1:3]<-lapply(plotS2[1:3],function(x) x+ylim(0,1.4))

  }else if(max(LMM_MSstat$Abundance,na.rm=T)<1.5){
    plotS2<-lapply(plotS2,function(x) x+ylim(0,1.25))

  }else if(max(target_protein1$Abundance,na.rm=T)<1.5){
    plotS2[1:3]<-lapply(plotS2[1:3],function(x) x+ylim(0,1.25))
  }

  plotS2<-lapply(plotS2,function(x) x+theme(text = element_text(size = 15),legend.position = "bottom",
                                            axis.title.x = element_text(size=20),
                                            axis.title.y=element_text(size=20),
                                            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=20),
                                            axis.text.y=element_text(size=20))+
                   scale_x_continuous("temperature",breaks=temps)+expand_limits(x = c(min(temps), max(temps))))
  P_fits<-ggpubr::ggarrange(plotlist=plotS2,ncol=4,nrow=1,font.label = list(size = 23),labels = l,legend.grob = y)
  return(P_fits)

}
