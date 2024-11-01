
Sim2MSstatsTMT<-function(Sim,half_temp_grid=FALSE){
  if(any(names(Sim)=="ProteinLevelData")){
    Sim<-Sim$ProteinLevelData |> as.data.frame()
  }
  if(isTRUE(half_temp_grid)){
    x<-Sim |> dplyr::mutate(
      BioReplicate=paste0(Protein,"_",temperature),
      Mixture=Protein,
      Run="Simulation",
      Experiment=Subject,
      Fraction=1,
      Condition=Group)#based on N=5 temps
  }else{
    x<-Sim |> dplyr::mutate(
      BioReplicate=paste0(Protein,"_",Channel),
      Mixture=Protein,
      Run="Simulation",
      Experiment=Subject,
      Fraction=1,
      Condition=Group)#based on N=10 temps

  }
  TechReps<-data.frame(Experiment=as.factor(unique(x$Subject)),TechRepMixture=rep(c(1,2),length(unique(x$Subject))/2))

  if(any(names(x)=="TechRepMixture")){
    x <-x |> dplyr::select(-TechRepMixture)
  }
  x<-x |> dplyr::inner_join(TechReps,by="Experiment")|> as.data.frame()
  return(x)
}
