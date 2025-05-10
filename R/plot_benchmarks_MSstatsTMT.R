plot_benchmarks_MSstatsTMT<-function(result,design="TPP",shifter="Non"){
  if (!requireNamespace("MSstatsTMT", quietly = TRUE)) {
    stop("The MSstatsTMT package is required but not installed.")
  }
  if (!requireNamespace("MSstats", quietly = TRUE)) {
    stop("The MSstats package is required but not installed.")
  }
  if (!requireNamespace("MSstatsConvert", quietly = TRUE)) {
    stop("The MSstatsConvert package is required but not installed.")
  }

  #QC plot of the simulation with 5 icc values
  #set temperature
  temps<-as.numeric(unique(result$temperature))
  #define an icc column based on the protein ID
  result$ICC<-stringr::str_extract(result$Protein,"icc_[:digit:].[[:digit:]]+")

  if(any(result$Condition==1)){
    result$treatment<-ifelse(result$Condition==1,"vehicle","treated")
    result$Condition<-result$treatment
    result$temperature<-53
  }
  if(design=="onePot"){
    result$Mixture<-"OnePot"
    result$Run<-"OnePot"
    result$Condition<-stringr::str_extract(result$Condition,"[[:lower:]]+")
  }

  annotation_file<-result|>dplyr::select(Run,Mixture,TechRepMixture,BioReplicate,Condition,Subject)|>dplyr::distinct()
  write.csv(annotation_file,paste0("annotation_file",shifter,"_",design,".csv"))
  #Define a data frame with one protein sim per ICC value
  One_prot_ICC<-result|>
    dplyr::group_by(ICC)|>
    dplyr::filter(Protein %in% unique(result$Protein)[2])|>
    dplyr::mutate(Accession=Protein)|>
    dplyr::ungroup()

  png(filename = paste0("ProfilePlots_",shifter,design,"_MsstatsTMTproc.png"),
      width =1600, height = 1600, units = "px", pointsize = 12,
      res = 600,type ="cairo")
  Profile_plot<-ggplot(One_prot_ICC,mapping=aes(x=Condition,y=Abundance,color=treatment))+geom_point(size=5)+
    ylab(expression(log[2]~Abundance))+
    geom_step(size=1.1)+
    ggtitle(paste0("Simulation template: ",shifter, " interaction"))+
    scale_x_continuous("Temperature",labels=temps,breaks=temps)+
    facet_wrap(~c(ICC),nrow=1)+theme(text=element_text(size=8))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="bottom")
  Profile_plot
  dev.off()
  saveRDS(Profile_plot,paste0("Profile_plot_",shifter,design,"_MsstatsTMTproc.RDS"))
  ##annotation file

    annotation_file<-result|>dplyr::select(Run,Mixture,TechRepMixture,BioReplicate,Condition,Subject,Channel)|>dplyr::distinct()
    annotation_file$Channel<-set_temps(10,temperatures=seq(1,10))$Channel

  annotation_file<-result|>dplyr::select(Run,Mixture,TechRepMixture,BioReplicate,Condition,Subject,Channel)|>dplyr::distinct()
  write.csv(annotation_file,paste0("annotation_file",shifter,"_",design,".csv"))

  dataMSstat<-list(ProteinLevelData=result)

  comparison<-make_contrast_matrix_all(dataMSstat,temps=NA)

  ATE_MSstats<-MSstatsTMT::groupComparisonTMT(
    dataMSstat,
    contrast.matrix = comparison,
    moderated = FALSE,
    adj.method = "BH",
    remove_norm_channel = TRUE,
    remove_empty_channel = FALSE,
    save_fitted_models = TRUE,
    use_log_file = FALSE,
    append = FALSE,
    verbose = TRUE,
    log_file_path = NULL
  )
  ATE_MSstats$ComparisonResult$ICC<-stringr::str_extract(ATE_MSstats$ComparisonResult$Protein,"icc_[:digit:].[[:digit:]]+")
  ATE_MSstats$ComparisonResult$ICC<-paste0("% of bio var = ",100*as.numeric(stringr::str_extract(ATE_MSstats$ComparisonResult$ICC,"[:digit:].[[:digit:]]+")))
  ATE_MSstats$ComparisonResult<-ATE_MSstats$ComparisonResult[stringr::str_detect(ATE_MSstats$ComparisonResult$ICC,"5|40"),]
  ATE_MSstats$ComparisonResult<-ATE_MSstats$ComparisonResult|>dplyr::group_by(ICC)|>dplyr::mutate(Sens=100*sum(pvalue<0.001)/length(unique(Protein)))
  png(filename = paste0("Histogram_",shifter,design,"_MsstatsTMTproc.png"),
      width =600, height = 600, units = "px", pointsize = 12,
      res = 130,type ="cairo")
  if(design=="onePot"){
    if(shifter=="strong"){
    MSstat_hist<-ggplot2::ggplot(ATE_MSstats$ComparisonResult,mapping=aes(x=pvalue))+
      geom_histogram(fill="#030366",color="black",bins=1,binwidth = 0.025)+facet_wrap(~factor(ICC,levels=c("% of bio var = 5","% of bio var = 40")),nrow=1)+
      coord_cartesian(xlim = c(0, 1))+ylim(0,1000)+xlab("pvalue")+
      theme(text=element_text(size=15),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+coord_cartesian(xlim = c(0, 1))+
      scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
      geom_text(mapping=aes(x=0.5,y=800),
                label=paste0(ATE_MSstats$ComparisonResult$Sens," %"),size=6)+ylab("protein count")
    }else{
      MSstat_hist<-ggplot2::ggplot(ATE_MSstats$ComparisonResult,mapping=aes(x=pvalue))+
        geom_histogram(fill="#030366",color="black")+facet_wrap(~factor(ICC,levels=c("% of bio var = 5","% of bio var = 40")),nrow=1)+
        scale_x_continuous(n.breaks=8)+ylim(0,1000)+xlab("pvalue")+
        theme(text=element_text(size=15),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+coord_cartesian(xlim = c(0, 1))+
        scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
        geom_text(mapping=aes(x=0.5,y=800),
                  label=paste0(ATE_MSstats$ComparisonResult$Sens," %"),size=6)+ylab("protein count")
    }
  }else{#if this is a TPP design
  MSstat_hist<-ggplot2::ggplot(ATE_MSstats$ComparisonResult,mapping=aes(x=pvalue))+
    geom_histogram(fill="#2C7FB8",color="black")+facet_wrap(~factor(ICC,levels=c("% of bio var = 5","% of bio var = 40")),nrow=1)+
    ylim(0,1000)+xlab("pvalue")+ scale_x_continuous(n.breaks=8)+
    theme(text=element_text(size=15),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+coord_cartesian(xlim = c(0, 1))+
    scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
    geom_text(mapping=aes(x=0.5,y=800),
              label=paste0(ATE_MSstats$ComparisonResult$Sens," %"),size=6)+ylab("protein count")

  }
  MSstat_hist
  dev.off()
  saveRDS(MSstat_hist,paste0("Histogram_",shifter,design,"_MsstatsTMTproc.RDS"))
  return(ATE_MSstats)
}
