plot_benchmarks_MSstatsTMT_reps<-function(result,design="TPP",shifter="Non",n_replicates_per_plex=10){
  #QC plot of the simulation with 5 icc values
  #define an icc column based on the protein ID
  result$ICC<-stringr::str_extract(result$Protein,"icc_[:digit:].[[:digit:]]+")
  if(any(result$Condition==1)){
    result$treatment<-ifelse(result$Condition==1,"vehicle","treated")

    result$Condition<-result$treatment
    result$temperature<-53
    if(length(unique(result$Subject))==20){
      result$Mixture<-result$Condition
    }
  }else{
    result$treatment<-result$Condition
  }
  if(n_replicates_per_plex==10){
    result<-result|>
      dplyr::group_by(Protein,Condition,temperature)|>
      dplyr::mutate(treatment=Condition)|>
      dplyr::group_split()|>
      lapply(function(x)x|>dplyr::mutate(replicate=seq(1,10)))|>
      dplyr::bind_rows()
  }
  if(design=="onePot"){
    if(n_replicates_per_plex==2|n_replicates_per_plex==5){
      result$Mixture<-"onePot"
      result$Run<-"onePot"
      result$Condition<-stringr::str_extract(result$Condition,"[[:lower:]]+")

      result<-result|>dplyr::group_by(Protein,Condition,replicate)|>dplyr::group_split()
    }else{
      result$Mixture<-stringr::str_extract(result$Condition,"[[:lower:]]+")
      result$Run<-result$Mixture
      result<-result|>dplyr::group_by(Protein,treatment,replicate)|>dplyr::group_split()

    }

    #unlog average and log back
    result<-lapply(result,function(x) {
      x$Abundance<-2^(x$Abundance)
      x$Abundance<-mean(x$Abundance,na.rm=T)
      x$Abundance<-log2(x$Abundance)
      if(design=="onePot"&n_replicates_per_plex==2|n_replicates_per_plex==5){
       x$BioReplicate=x$BioReplicate[1]
       x$Subject=x$Subject[1]
       x$Channel=x$Channel[1]
       x<-x|>dplyr::distinct()
      }
      return(x)
    })|>dplyr::bind_rows()
  }
  if(n_replicates_per_plex==5&design=="onePot"){
    annotation_file<-result|>dplyr::select(Run,Mixture,TechRepMixture,BioReplicate,Condition,Subject,Channel)|>dplyr::distinct()
    annotation_file$Channel<-set_temps(10,temperatures=seq(1,10))$Channel
    result<-result|>dplyr::select(-Channel)|>dplyr::inner_join(annotation_file)
  }else if(n_replicates_per_plex==2&design=="onePot"){
    annotation_file<-result|>dplyr::select(Condition,Subject,temperature)|>dplyr::distinct()
    annotation_file$Channel<-set_temps(10,temperatures=seq(1,10))$Channel[1:4]
    result<-result|>dplyr::select(-Channel)|>dplyr::inner_join(annotation_file)
  }
  annotation_file<-result|>dplyr::select(Run,Mixture,TechRepMixture,BioReplicate,Condition,Subject,Channel)|>dplyr::distinct()
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
    geom_step(size=1.1)+
    ylab(expression(log[2]~Abundance))+
    ggtitle(paste0("Simulation template: ",shifter, " interaction"))+
    scale_x_continuous("Temperature",breaks=as.numeric(unique(temps)), labels=as.numeric(unique(temps)))+
    facet_wrap(~c(ICC),nrow=1)+theme(text=element_text(size=15))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="bottom")
  Profile_plot
  dev.off()
  saveRDS(Profile_plot,paste0("Profile_plot_",shifter,design,"_MsstatsTMTproc.RDS"))
  dataMSstat<-list(ProteinLevelData=result)

  comparison<-make_contrast_matrix_all(dataMSstat,temps=NA)

  ATE_MSstats<-MSstatsTMT::groupComparisonTMT(
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
  return(ATE_MSstats)
}
