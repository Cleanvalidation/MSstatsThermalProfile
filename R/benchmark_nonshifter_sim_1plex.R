
benchmark_nonshifter_sim_1plex<-function(msstats_icc_output,templateProtein,n_sims,t_range=seq(7,7),design="OnePot",n_replicates_per_plex=NA){
  if (!requireNamespace("MSstatsTMT", quietly = TRUE)) {
    stop("The MSstatsTMT package is required but not installed.")
  }
  if (!requireNamespace("MSstats", quietly = TRUE)) {
    stop("The MSstats package is required but not installed.")
  }
  if (!requireNamespace("MSstatsConvert", quietly = TRUE)) {
    stop("The MSstatsConvert package is required but not installed.")
  }
  temps<-as.character(sort(unique(msstats_icc_output$df_with_variance$temperature))[t_range])
  #all proteins is the normalized processed output from msstats_icc
  #msstats_icc_output<-msstats_icc(MSstats_Humanproc_wImputation,temps=unique(MSstats_Humanproc_wImputation$ProteinLevelData$temperature))
  all_proteins<-msstats_icc_output$df_with_variance|>
    dplyr::mutate(Protein = as.character(Protein),
                  temperature= as.character(temperature))

  template_MSstats<-list(ProteinLevelData=all_proteins|>
                           dplyr::filter(Protein %in% templateProtein,
                                         temperature %in% temps))

  #get predicted values per condition

  #rename columns of the template to get population-level abundance values

  template_MSstats$ProteinLevelData$LogIntensities = template_MSstats$ProteinLevelData$Abundance

  template_MSstats$ProteinLevelData$GROUP = template_MSstats$ProteinLevelData$Condition


  template_MSstats$ProteinLevelData$SUBJECT = template_MSstats$ProteinLevelData$BioReplicate


  template_MSstats<-MSstats::quantification(template_MSstats, type="Group", use_log_file = FALSE)


  template_MSstats<-template_MSstats|>tidyr::pivot_longer(cols=colnames(template_MSstats)[colnames(template_MSstats)!="Protein"],
                                                          names_to="Condition",
                                                          values_to="Abundance")



  #set th template for hte simulation
  template_simulation<-template_MSstats
  #define mapping between temperatures and TMT channels
  temps<-all_proteins|>
    dplyr::mutate(Protein %in% templateProtein)|>
    dplyr::select(Channel,temperature)|>
    dplyr::distinct()

  #select a template protein

  template_simulation<-template_simulation|>
     tidyr::separate(Condition,into=c("Channel","treatment"))|>
    dplyr::mutate(Condition=paste0(Channel,"_",treatment),
                  BioReplicate=Condition)|>
    dplyr::inner_join(temps)|>dplyr::group_by(temperature)|>
    dplyr::mutate(Abundance=mean(Abundance,na.rm=T))

  temps<-sort(unique(all_proteins$temperature))[t_range]

  png(filename = "template_MsstatsTMTproc_nonshifter.png",
      width =12, height = 12, units = "in", pointsize = 12,
      res = 600,type ="cairo")
  Template<-ggplot(template_simulation,mapping=aes(x=treatment,y=Abundance,color=treatment))+
    geom_point(size=4)+
    # geom_step(linewidth=1)+
    # scale_x_continuous("Temperature", breaks = unique(template_simulation$temperature))+
    labs(title="Simulation template: non interaction",x="Temperature",y="Log of protein abundances")+theme_bw()+
    theme(text=element_text(size=20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  Template
  dev.off()
  saveRDS(Template,"templateMSstatsTMTproc_nonshifter.RDS")
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

  if(n_replicates_per_plex ==5){
    #number of replicates per condition = 5 because it is 1-plex
    n_replicates<-5
  }else if (n_replicates_per_plex == 2){
    n_replicates<-2
  }

  df<-template_simulation|>
    dplyr::select(Abundance,temperature,treatment,Channel)|>
    dplyr::arrange(temperature)
  df$Channel<- Channels[seq(1:nrow(df))]
  result<-list()
  #qc template
  #ggplot(template,mapping=aes(x=temperature,y=Abundance,color=treatment))+geom_point()
  for(i in seq(n_sims)){#n_sims is the number of simulations defined by the user

    for(j in icc_range){#for each protein, one intra-class correlation is used to generate a profile

      #Simulate profiles with variation components
      Sims<- add_variation_shift_1plex(df,j,n_conditions,n_replicates,re_var,t_range,design=design)

      #Revert the list back to a data frame
      Sims<-dplyr::bind_rows(Sims)

      #combine the simulation number with the ICC value
      Sims$Protein<-paste0("Sim_",i,"_icc_",j)

      #append simulated profile with the ICC to the output
      result<- append(result,list(Sims))

    }

  }

  #convert list data frames for sims per ICC value to data frame
  clean_result <- result[sapply(result, inherits, what = "data.frame")]
  result <- dplyr::bind_rows(clean_result)
  if(design=="OnePot"){
    if(n_replicates_per_plex==5){
      result<-result|>
        dplyr::group_by(Protein,Condition)|>
        dplyr::mutate(Channel = unique(msstats_icc_output$df_with_variance$Channel)[Subject],
                      temperature = mean(as.numeric(template_simulation$temperature),na.rm=T),
                      Run = "OnePot",
                      Mixture = "OnePot")|>
        dplyr::group_split()|>
        purrr::map(function(x) x[1:5,])|>
        dplyr::bind_rows()|>
        dplyr::ungroup()
    }else if (n_replicates_per_plex==2){
      result<-result|>
        dplyr::group_by(Protein,Condition)|>
        dplyr::mutate(Channel = unique(msstats_icc_output$df_with_variance$Channel)[Subject],
                      temperature = mean(as.numeric(template_simulation$temperature),na.rm=T),
                      Run = "OnePot",
                      Mixture = "OnePot")|>
        dplyr::group_split()|>
        purrr::map(function(x) x[1:2,])|>
        dplyr::bind_rows()|>

        dplyr::ungroup()
    }
  }

  return(result)
  #output sigmoid simulations
}

