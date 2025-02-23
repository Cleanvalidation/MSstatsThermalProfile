plot_benchmarksTPP<-function(result,temps=unique(result$temperature),shifter="strong",design="TPP"){
  #QC plot of the simulation with 5 icc values
  #define an icc column based on the protein ID
  result$ICC<-stringr::str_extract(result$Protein,"icc_[:digit:].[[:digit:]]+")
  if(design=="TPP"){
  result$Mixture<-result$Subject
  result$Run<-result$Subject
  }else{
    result$Run<-1
    result$Mixture<-"OnePot"
    result$TechRepMixture<-1

  }
  if(any(result$Condition==1)){
    result$treatment<-ifelse(result$Condition==1,"vehicle","treated")
    result$Condition<-result$treatment
    result$temperature<-53
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
      res = 130,type ="cairo")
  Profile_plot<-ggplot(One_prot_ICC,mapping=aes(x=Condition,y=Abundance,color=treatment))+geom_point()+
    ylab(expression(log[2]~Abundance))+
    facet_wrap(~c(ICC),nrow=1)+theme(text=element_text(size=15))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  Profile_plot
  dev.off()
  saveRDS(Profile_plot,paste0("Profile_plot_",shifter,design,"_MsstatsTMTproc.RDS"))
  dataTPP<-list(ProteinLevelData=result)
  comparison<-make_contrast_matrix(dataMSstat,temps=temps)
  #set condition vehicle and treated for TPP
 dataTPP$ProteinLevelData$Condition<-stringr::str_extract(dataTPP$ProteinLevelData$treatment,"[[:lower:]]+")
 #set technical replicate
 dataTPP$ProteinLevelData$TechRepMixture<-stringr::str_extract(dataTPP$ProteinLevelData$Subject,"_[:digit:]")
 #Run data with TPP splines
  TPP<-TPP_NPARC_calc(dataTPP$ProteinLevelData,method="NPARC",DF=5,CARRIER=FALSE,temps=set_temps(10,c(37.3, 40.6, 43.9, 47.2, 50.5, 53.8, 57.1, 60.4, 64, 67)),NORM=FALSE,filters=TRUE)


    MSstatsTMT::groupComparisonTMT(
    dataMSstat,
    contrast.matrix = comparison,
    moderated = FALSE,
    adj.method = "BH",
    remove_norm_channel = TRUE,
    remove_empty_channel = TRUE,
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

}
