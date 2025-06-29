#' @importFrom stringr str_extract
plot_benchmarks_MSstatsTMT_reps<-function(result,design="TPP",shifter="Non",t_range=seq(1,10), n_replicates_per_plex=10,variation_idx=NA){
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
  #define an icc column based on the protein ID
  result$ICC<-stringr::str_extract(result$Protein,"icc_[:digit:].[[:digit:]]+")
  if(any(result$Condition==1)){#If this is a OnePot design, assign one temperature by default
    result$treatment<-ifelse(result$Condition==1,"vehicle","treated")

    result$Condition<-result$treatment
    result$temperature<-53
    if(length(unique(result$Subject))==20){
      result$Mixture<-result$Condition
    }
  }else{
    result$treatment<-result$Condition
  }
  if(n_replicates_per_plex==4&length(t_range)==2){
    temps<-unique(result$temperature)[t_range] #the temperatures for the contrast
    #only keep the temperatures for the contrast
    result<-result|>
      dplyr::filter(temperature %in% temps)
  }
  if(n_replicates_per_plex==10){
    result<-result|>
      dplyr::group_by(Protein,Condition,temperature)|>
      dplyr::mutate(treatment=Condition)|>
      dplyr::group_split()|>
      lapply(function(x)x|>dplyr::mutate(replicate=seq(1,10)))|>
      dplyr::bind_rows()
  }
  temps<-unique(result$temperature)
  if(design=="OnePot"){
    if(n_replicates_per_plex==5){
      annotation_file<-result|>
        dplyr::select(Run,Mixture,TechRepMixture,BioReplicate,Condition,Subject,temperature, treatment, Channel)|>
        dplyr::distinct()|>
        dplyr::group_by(Condition)|>
        dplyr::mutate(Replicate=seq(1,10),
                      Subject=paste0(ifelse(Condition=="vehicle",1,2),"_",Replicate),
                      BioReplicate=Subject,
                      Condition=treatment,
                      Run = paste0("OnePot",sample(1:2,10,replace=TRUE)),
                      Mixture = Run)|>
        dplyr::ungroup()|>
        dplyr::select(-Replicate)

      result<-result|>
        dplyr::select(-Subject,-BioReplicate,-Mixture,-Run)|>
        dplyr::inner_join(annotation_file)

    }
  }
 # Keep 2 replicates per protein at random
 if(n_replicates_per_plex==2|n_replicates_per_plex==5){
   result<-result|>
     dplyr::group_by(Protein,Condition)|>
     dplyr::group_split()|>
     lapply(function(x) x[sample(1:nrow(x),size=n_replicates_per_plex),])|>
     dplyr::bind_rows()|>
     dplyr::ungroup()
 }else{

  annotation_file<-result|>
    dplyr::select(Run,Mixture,TechRepMixture,BioReplicate,Condition,Subject,Channel)|>dplyr::distinct()
  write.csv(annotation_file,paste0("annotation_file",shifter,"_",design,".csv"))
 }
  #Define a data frame with one protein sim per ICC value
  One_prot_ICC<-result|>
    dplyr::group_by(ICC)|>
    dplyr::filter(Protein %in% unique(result$Protein)[2])|>
    dplyr::mutate(Accession=Protein)|>
    dplyr::ungroup()
  png(filename = paste0("ProfilePlots_",shifter,design,"_MsstatsTMTproc.png"),
      width =12, height = 6, units = "in", pointsize = 12,
      res = 600,type ="cairo")
  Profile_plot<-ggplot(One_prot_ICC,mapping=aes(x=temperature,y=Abundance,color=treatment))+geom_point()+
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

  if(length(t_range)==2&&any(!is.na(variation_idx))){
    temps<-unique(dataMSstat$ProteinLevelData$temperature)[t_range]
    variation_temps<-temps
    dataMSstat$ProteinLevelData <- dataMSstat$ProteinLevelData[dataMSstat$ProteinLevelData$temperature %in% variation_temps,]

  }else if(length(t_range)==2&&any(is.na(variation_temps))){
    temps<-unique(dataMSstat$ProteinLevelData$temperature)[t_range]
  }else{
    temps <- NA

  }
  comparison<-make_contrast_matrix_all(dataMSstat,temps=temps,variation_temps=variation_temps)


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
      width =12, height = 6, units = "in", pointsize = 12,
      res = 600,type ="cairo")
  if(design=="OnePot"){
    if(shifter=="strong"){
    MSstat_hist<-ggplot2::ggplot(ATE_MSstats$ComparisonResult,mapping=aes(x=pvalue))+
      geom_histogram(fill="#030366",color="black",bins=1,binwidth = 0.025)+facet_wrap(~factor(ICC,levels=c("% of bio var = 5","% of bio var = 40")),nrow=1)+
      coord_cartesian(xlim = c(0, 1))+ylim(0,1000)+xlab("pvalue")+
      theme(text=element_text(size=12),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      coord_cartesian(xlim = c(0, 1))+
      scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
      geom_text(mapping=aes(x=0.5,y=800),
                label=paste0(ATE_MSstats$ComparisonResult$Sens," %"),size=8)+ylab("protein count")
    }else{
      MSstat_hist<-ggplot2::ggplot(ATE_MSstats$ComparisonResult,mapping=aes(x=pvalue))+
        geom_histogram(fill="#030366",color="black")+facet_wrap(~factor(ICC,levels=c("% of bio var = 5","% of bio var = 40")),nrow=1)+
        scale_x_continuous(n.breaks=8)+ylim(0,1000)+xlab("pvalue")+
        theme(text=element_text(size=12),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+coord_cartesian(xlim = c(0, 1))+
        scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
        geom_text(mapping=aes(x=0.5,y=800),
                  label=paste0(ATE_MSstats$ComparisonResult$Sens," %"),size=8)+ylab("protein count")
    }
  }else{#if this is a TPP design
  MSstat_hist<-ggplot2::ggplot(ATE_MSstats$ComparisonResult,mapping=aes(x=pvalue))+
    geom_histogram(fill="#2C7FB8",color="black")+facet_wrap(~factor(ICC,levels=c("% of bio var = 5","% of bio var = 40")),nrow=1)+
    ylim(0,1000)+xlab("pvalue")+ scale_x_continuous(n.breaks=8)+
    theme(text=element_text(size=12),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+coord_cartesian(xlim = c(0, 1))+
    scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
    geom_text(mapping=aes(x=0.5,y=800),
              label=paste0(ATE_MSstats$ComparisonResult$Sens," %"),size=6)+ylab("protein count")

  }
  MSstat_hist
  dev.off()
  saveRDS(MSstat_hist,paste0("Histogram_",shifter,design,"_MsstatsTMTproc.RDS"))
  return(ATE_MSstats)
}
