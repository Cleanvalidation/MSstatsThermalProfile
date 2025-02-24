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
    dplyr::mutate(Accession=Protein,
                  Condition=stringr::str_extract(treatment,"[[:lower:]]+"))|>
    dplyr::ungroup()
  png(filename = paste0("ProfilePlots_",shifter,design,"_TPPproc.png"),
      width =1600, height = 1600, units = "px", pointsize = 12,
      res = 600,type ="cairo")
  Profile_plot<-ggplot(One_prot_ICC,mapping=aes(x=temperature,y=Abundance,color=treatment))+
    geom_point()+geom_step(size=1.1)+
    ylab("Log of protein abundances")+
    ggtitle(paste0("Simulation template: ",shifter, " interaction"))+
    scale_x_continuous("Temperature",breaks=as.numeric(unique(result$temperature)), labels=as.numeric(unique(result$temperature)))+
    facet_wrap(~c(ICC),nrow=1)+theme(text=element_text(size=8))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  Profile_plot
  dev.off()
  saveRDS(Profile_plot,paste0("Profile_plot_",shifter,design,"_TPPproc.RDS"))
  dataTPP<-list(ProteinLevelData=result)
  comparison<-make_contrast_matrix(dataMSstat,temps=temps)
  #set condition vehicle and treated for TPP
 dataTPP$ProteinLevelData$Condition<-stringr::str_extract(dataTPP$ProteinLevelData$treatment,"[[:lower:]]+")
 #set technical replicate
 dataTPP$ProteinLevelData$TechRepMixture<-stringr::str_extract(dataTPP$ProteinLevelData$Subject,"_[:digit:]")
 #reference channel
 refChannel<-dataTPP$ProteinLevelData$Channel[which(unique(dataTPP$ProteinLevelData$temperature)==min(unique(dataTPP$ProteinLevelData$temperature)))]
 #add necessary column names
 dataTPP$ProteinLevelData$Accession<-dataTPP$ProteinLevelData$Protein
 dataTPP$ProteinLevelData$File.ID<-dataTPP$ProteinLevelData$Subject
 #normalize with TPP dataTPP<-TPP_normalization(dataTPP$ProteinLevelData,TPPfilters=FALSE,temps=unique(temps),reference=refChannel,CARRIER=FALSE)
 #Run data with TPP splines
  TPP<-TPP_NPARC_calc(dataTPP$normData,method="NPARC",DF=5,CARRIER=FALSE,temps=set_temps(10,c(37.3, 40.6, 43.9, 47.2, 50.5, 53.8, 57.1, 60.4, 64, 67)),NORM=FALSE,filters=TRUE)

  TPP$ICC<-stringr::str_extract(TPP$uniqueID,"icc_[:digit:].[[:digit:]]+")
  #append biological variance text to numeric values
  TPP$ICC<-paste0("% of bio var = ",100*as.numeric(stringr::str_extract(TPP$ICC,"[:digit:].[[:digit:]]+")))
  #Keep 5 and 40% ICC
  TPP<-TPP[stringr::str_detect(TPP$ICC,"5|40"),]
  TPP<-TPP|>dplyr::group_by(ICC)|>dplyr::mutate(Sens=100*sum(p_adj_NPARC<0.001)/length(unique(TPP$uniqueID)))
  png(filename = paste0("Histogram_",shifter,design,"_TPPproc.png"),
      width =600, height = 600, units = "px", pointsize = 12,
      res = 130,type ="cairo")
  if(design=="onePot"){
    if(shifter=="strong"){
      TPP_hist<-ggplot2::ggplot(TPP,mapping=aes(x=p_adj_NPARC))+
        geom_histogram(fill="#D95F0E",color="black",bins=1,binwidth = 0.025)+facet_wrap(~factor(ICC,levels=c("% of bio var = 5","% of bio var = 40")),nrow=1)+
        coord_cartesian(xlim = c(0, 1))+ylim(0,1000)+xlab("pvalue")+
        theme(text=element_text(size=15),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+coord_cartesian(xlim = c(0, 1))+
        scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
        geom_text(mapping=aes(x=0.5,y=800),
                  label=paste0(TPP$Sens," %"),size=6)+ylab("protein count")
    }else{
      TPP_hist<-ggplot2::ggplot(TPP$ComparisonResult,mapping=aes(x=p_adj_NPARC))+
        geom_histogram(fill="#D95F0E",color="black")+facet_wrap(~factor(ICC,levels=c("% of bio var = 5","% of bio var = 40")),nrow=1)+
        scale_x_continuous(n.breaks=8)+ylim(0,1000)+xlab("pvalue")+
        theme(text=element_text(size=15),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+coord_cartesian(xlim = c(0, 1))+
        scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
        geom_text(mapping=aes(x=0.5,y=800),
                  label=paste0(TPP$Sens," %"),size=6)+ylab("protein count")
    }
  }else{#if this is a TPP design
    TPP_hist<-ggplot2::ggplot(TPP$ComparisonResult,mapping=aes(x=p_adj_NPARC))+
      geom_histogram(fill="#D95F0E",color="black")+facet_wrap(~factor(ICC,levels=c("% of bio var = 5","% of bio var = 40")),nrow=1)+
      ylim(0,1000)+xlab("pvalue")+ scale_x_continuous(n.breaks=8)+
      theme(text=element_text(size=15),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+coord_cartesian(xlim = c(0, 1))+
      scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0))+
      geom_text(mapping=aes(x=0.5,y=800),
                label=paste0(TPP$Sens," %"),size=6)+ylab("protein count")

  }
  TPP_hist
  dev.off()
  saveRDS(TPP_hist,paste0("Histogram_",shifter,design,"_TPPproc.RDS"))

}
