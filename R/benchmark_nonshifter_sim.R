benchmark_nonshifter_sim<-function(msstats_icc_output,templateProtein,n_sims,t_range=seq(7,7),design="TPP"){

  if (!requireNamespace("MSstatsTMT", quietly = TRUE)) {
    stop("The MSstatsTMT package is required but not installed.")
  }
  if (!requireNamespace("MSstats", quietly = TRUE)) {
    stop("The MSstats package is required but not installed.")
  }
  if (!requireNamespace("MSstatsConvert", quietly = TRUE)) {
    stop("The MSstatsConvert package is required but not installed.")
  }

  #all proteins is the normalized processed output from msstats_icc
  #msstats_icc_output<-msstats_icc(MSstats_Humanproc_wImputation,temps=unique(MSstats_Humanproc_wImputation$ProteinLevelData$temperature))
  all_proteins<-msstats_icc_output$df_with_variance
  template_MSstats<-list(ProteinLevelData=all_proteins|>
                           dplyr::filter(Protein %in% templateProtein))
  #get predicted values per condition
  #rename columns of the template to get population-level abundance values
  template_MSstats$ProteinLevelData$LogIntensities = template_MSstats$ProteinLevelData$Abundance
  template_MSstats$ProteinLevelData$GROUP = paste0(template_MSstats$ProteinLevelData$Channel,"_",template_MSstats$ProteinLevelData$treatment)
  template_MSstats$ProteinLevelData$SUBJECT = template_MSstats$ProteinLevelData$BioReplicate
  template_MSstats<-MSstats::quantification(template_MSstats, type="Group", use_log_file = FALSE)
  template_MSstats<-template_MSstats|>tidyr::pivot_longer(cols=colnames(template_MSstats)[colnames(template_MSstats)!="Protein"],
                                                          names_to="Condition",
                                                          values_to="Abundance")
  #Averaging vehicle and treated for each Channel so that abundance is equal for both the conditions
  template_MSstats$Channel<-stringr::str_split_fixed(template_MSstats$Condition,"_",2)[,1]
  data<-dplyr::select(template_MSstats,Channel,Abundance) %>%
    dplyr::group_by(Channel ) %>%
    tidyr::pivot_wider(names_from = Channel, values_from = Abundance,values_fn=mean)
  data<-data|>
    tidyr::pivot_longer(cols=colnames(data),
      names_to="Channel",values_to="Abundance")



  # Joining the averaged data with template_Msstats
  template_simulation<-template_MSstats|>dplyr::select(-Abundance)|>dplyr::inner_join(data)

  #define mapping between temperatures and TMT channels
  temps<-all_proteins|>dplyr::mutate(Protein %in% templateProtein)|>dplyr::select(Channel,temperature,treatment)|>dplyr::distinct()
  #select a template protein
  template_simulation<-template_simulation|>
    tidyr::separate(Condition,into=c("Channel","treatment"))|>
    dplyr::mutate(Condition=paste0(Channel,"_",treatment),
                  BioReplicate=Condition)|>
    dplyr::inner_join(temps)
  #save non-shifting template
  png(filename = "template_MsstatsTMTproc_nonshift.png",
      width =1600, height = 1600, units = "px", pointsize = 12,
      res = 600,type ="cairo")
  Template_TN<-ggplot(template_simulation,mapping=aes(x=temperature,y=Abundance,size=treatment,shape=treatment,color=treatment))+
    geom_point(size=4,alpha=0.5)+
    geom_step(size=1.1)+
    scale_size_manual(values=c(8,5))+
    ylab("Log of protein abundances")+
    ggtitle(paste0("Simulation template: non interaction"))+
    scale_x_continuous("Temperature",breaks=as.numeric(unique(result$temperature)), labels=as.numeric(unique(result$temperature)))+
    theme(text=element_text(size=8))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="bottom")
  Template_TN
  dev.off()
  saveRDS(Template_TN,"Template_TN.RDS")

  #define inputs
  # dplyr::group_split()|>
  # dplyr::first()
  icc_range<-c(0.05,0.2,0.4,0.6,0.8)#from the histogram data, select a range of icc values (0-0.8)
  re_var<-median(all_proteins$sigma_re^2)
  #define output
  result<-list()
  #define the number of conditions
  n_conditions<-2
  #define the number of temps
  n_temps<-length(t_range)
  #define the number of replicates per condition based on the number of temperatures provided
  #number of replicates per condition
  if(n_temps==1){
    n_replicates<-20
    #template_simulation$Subject<-paste0(template_simulation$temperature,"_",template_simulation$TechRepMixture,"_",paste0(template_simulation$Condition))#organize this in one place
  }else if(n_temps==2){
    n_replicates<-10
    #template_simulation$Subject<-paste0(template_simulation$temperature,"_",template_simulation$TechRepMixture,"_",paste0(template_simulation$Condition))#organize this in one place
  }else if(n_temps==5){
    n_replicates<-4
    #template_simulation$Subject<-paste0(template_simulation$temperature,"_",template_simulation$TechRepMixture,"_",paste0(template_simulation$Condition))#organize this in one place
  }else if(n_temps==10){
    n_replicates<-2
  }
  df<-template_simulation|>
    dplyr::select(Abundance,temperature,treatment,Channel)
  result<-list()
  #qc template
  #ggplot(template,mapping=aes(x=temperature,y=Abundance,color=treatment))+geom_point()
  for(i in seq(n_sims)){#n_sims is the number of simulations defined by the user
    for(j in icc_range){#for each protein, one intra-class correlation is used to generate a profile
      #Simulate profiles with variation components
      Sims<- add_variation_shift(df,j,n_conditions,n_replicates,re_var,t_range,design=design)
      #Revert the list back to a data frame
      Sims<-dplyr::bind_rows(Sims)
      # #Map the Sims subject, condition and treatment indices to the actual values in the template data
      #       #if this is a 1D-TPP assay, define Condition columns
      #combine the simulation number with the ICC value
      Sims$Protein<-paste0("Sim_",i,"_icc_",j)
      #append simulated profile with the ICC to the output
      result<- append(result,list(Sims))
    }
  }
  #convert list data frames for sims per ICC value to data frame
  result<-dplyr::bind_rows(result)
  return(result)
  #output sigmoid simulations
}
