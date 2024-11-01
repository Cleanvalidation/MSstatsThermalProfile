add_variation_shift<-function(template,icc,n_conditions,n_replicates,re_var,t_range,design="TPP"){
  n_temps<-length(t_range)
  #Define subjects per condition
  if(n_temps==1){
    n_subjects<-40
  }else if(n_temps==2){
    n_subjects<-10
  }else if(n_temps==5){
    n_subjects<-4
  }else if(n_temps==10){
    n_subjects<-2
  }else{#if the number of temps is not within specified range, stop the script
    stopifnot("Please define n_temps as 1, 2, 5 or 10" = n_temps %in% c(1,2,5,10))
  }

  #QC template curve
  #ggplot(template,mapping=aes(y=Abundance,x=temperature,color=treatment))+geom_point()
  #df is a vector of n_temps abundances corresponding to n_temps (input is per condition)
  ##for shifters df is a matrix 10x2 the first column is the abundances and the second column is the condition
  template$treatment<-tolower(template$treatment)
  template<-template|>dplyr::mutate(Condition=ifelse(treatment=="vehicle",
                                                     as.character(1),
                                                     as.character(2)))

  # template_condition<-template|>
  #   dplyr::select(Abundance,Condition,temperature)|>dplyr::distinct()
  #
  # #define input as a 10x2
  # input<-template_condition|>tidyr::pivot_wider(names_from="Condition",values_from=Abundance)|>dplyr::select(-temperature)
  #

  #first column prediction for vehicle
  #second column prediction for treated
  #icc is a constant intra-class correlation
  #  #re_var is a constant residual error variance
  #  #n_protein is the number of proteins simulated
  #n_conditions
  #n_replicates

  n_obs<-n_temps*n_conditions*n_replicates
  stopifnot("number of observations must be 40 for the experimental design with 4 TMT-10plexes" = n_obs==40)

  #modified_matrix[,3]<-seq(1:n_temps)

  biovar<-(re_var)/(1/icc-1)

  #select  the abundances at the specified temperatures and replicate the fitted values for each TMT 10-plex and each condition
  if(n_temps==1 & design!="hybrid"){#in this case the input has 1 row and 40 cols, each column being a separate subject
    #take the temperature index from the input
    template<-template|>dplyr::group_by(Condition,TechRepMixture)|>dplyr::group_split()
    template<-lapply(template,function(x) x|>dplyr::mutate(Abundance=rep(x$Abundance[t_range],n_replicates/n_conditions),
                                                           temperature=rep(x$temperature[t_range],n_replicates/n_conditions)))|>
      dplyr::bind_rows()


    # create an empty matrix to store modified responses
    modified_matrix <- matrix(nrow = 1, ncol = 40)

    for ( i in 1:n_conditions){
      #iterate over the number of subjects defined as an input and generate biological variation per subject
      for (j in 1:n_subjects){
        sample_u<-rnorm(1,mean=0,sd=sqrt(biovar))
        current_responses<-template$Abundance[template$Subject==unique(template$Subject)[j]]
        #iterate over the temperatures selected
        for (k in 1:n_temps){
          #sample residual error per temperature
          sample_e<-rnorm(1,mean=0,sd=sqrt(re_var))
          #get Y_true
          Y_true<- as.numeric(current_responses)

          Y_observed<- Y_true+sample_u+sample_e

          # # store modified responses in new matrix
          modified_matrix[k,j]<-Y_observed#10x4 matrix for TMT 10-plex, each column indicating the subject

        }

      }
    }


    colnames(modified_matrix)<-paste0(template$temperature,"_",rownames(template),"_",template$treatment)
    modified_df<-as.data.frame(modified_matrix)

    #
    result<-tidyr::pivot_longer(modified_df,
                                cols=colnames(modified_df)[!colnames(modified_df)=="Temperature"],
                                names_to="Subject",
                                values_to="Abundance")|>
      dplyr::mutate(Condition=ifelse(stringr::str_detect(Subject,"vehicle"),1,2),
                    TechRepMixture=1)
    #generate temperature values selected


    result$Temperature<-rep(t_range,nrow(result))

    result$BioReplicate<-result$Subject
    result$Mixture<-result$Subject
    result$Run<-result$Subject
    result$treatment<-ifelse(result$Condition=="1","vehicle","treatment")
    result$temperature<-stringr::str_remove(stringr::str_extract(result$Subject,"[[:digit:]]+.[:digit:]|[[:digit:]]+_"),"_")
    result$Condition<-paste0(result$temperature,"_",result$treatment)
    result$Channel<-unique(template$Channel)[t_range]
    result$TechRepMixture<-1
  }else if(n_temps==1 & design=="hybrid"){#If this is a one-pot hybrid design
    n_temps<-10
    n_subjects<-1
    # create an empty matrix to store modified responses

    one_pot_list<-list()
    #modified matrix is 1 row,and 2 columns, one column per condition
    modified_matrix <- matrix(NA, nrow = 10, ncol = n_conditions)

    temperatures<-unique(template$temperature)[t_range]
    #input is a 10x2 data frame
    template<-template|>dplyr::select(Condition,Abundance,Channel)|>dplyr::distinct()

    for (i in 1:n_conditions) {
      for (j in 1:n_subjects) {

        sample_u <- rnorm(1, mean = 0, sd = sqrt(biovar))#simulate subject effect

        for (k in 1:n_temps) {
          current_responses <- template$Abundance[template$Condition == i][k]
          sample_e <- rnorm(1, mean = 0, sd = sqrt(re_var))#simulate random error
          Y_true <- as.numeric(current_responses)
          Y_observed <- Y_true + sample_u + sample_e#add biological variance and error term
          #for a separate function, unlog Y_observed, average Y_obs for all temps, log
          modified_matrix[k , i] <- Y_observed
        }
      }

      #data frame that averages abundances from all channels per condition
      one_pot_df<-data.frame(Condition=i,Abundance=mean(modified_matrix[,i],na.rm=T))
      #append results to a list, each index corresponding to one condition
      one_pot_list<-append(one_pot_list,list(one_pot_df))
    }

    #append all conditions to one data frame
    result<-dplyr::bind_rows(one_pot_list)


    #generate temperature values selected
    result$Temperature<-rep(t_range,nrow(result))
    result$Subject<-paste0("one_pot_",unique(temperatures),"_",result$Condition)
    result$BioReplicate<-result$Subject
    result$Mixture<-result$Subject
    result$Run<-result$Subject
    result$treatment<-ifelse(result$Condition=="1","vehicle","treatment")
    result$temperature<-stringr::str_remove(stringr::str_extract(result$Subject,"[[:digit:]]+.[:digit:]|[[:digit:]]+_"),"_")
    result$Condition<-paste0(result$temperature,"_",result$treatment)
    result$Channel<-unique(template$Channel)[5]
    result$TechRepMixture<-1
  }else if(n_temps>1){# if the design has more than one temperature, add variance components

    result<-add_variance_components(template,biovar,re_var,n_conditions,n_subjects,n_temps,t_range)

  }

  return(result)
}
