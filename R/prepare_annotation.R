prepare_annotation<-function(template,modified_df){

  modified_df <- tidyr::pivot_longer(modified_df,
                                     cols = names(modified_df)[!names(modified_df) %in% c("temperature","Subject","Condition")],
                                     names_to = "Condition",
                                     values_to = "Abundance")



  modified_df$TechRepMixture<-1

  modified_df$treatment<-ifelse(modified_df$Condition=="1","vehicle","treatment")

  modified_df$Subject<-paste0(modified_df$Condition,"_",modified_df$Subject)

  modified_df$Condition<-paste0(modified_df$temperature,"_",modified_df$treatment)

  modified_df$Run<-modified_df$Subject

  modified_df$Mixture<-modified_df$Subject

  modified_df$BioReplicate<-modified_df$Subject
  if(any(names(template)=="Channel")){
    df<-subset(template, select = c("Channel","temperature"))
    df<-as.data.frame(unique(df))
    modified_df<-merge(modified_df, df, by='temperature', all.x=TRUE)
    annotation_file<-modified_df|>dplyr::select(Run,TechRepMixture,Channel,Condition,Mixture)


  }else{
    annotation_file<-data.frame()
  }
  #return a list result with the modified data frame format and annotation file
  result<-list(Simulation=modified_df,annotation_file=annotation_file)

  return(result)
}

