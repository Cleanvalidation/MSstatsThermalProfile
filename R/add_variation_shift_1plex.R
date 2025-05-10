#' @importFrom stringr str_extract


add_variation_shift_1plex<-function(template,icc,n_conditions,n_replicates,re_var,t_range,design="TPP"){
  n_temps<-length(t_range)
  #Define subjects per condition
  if(n_temps==1){
    n_subjects<-10
  }else if(n_temps==2){
    n_subjects<-5
  }else if(n_temps==5){
    n_subjects<-2
  }else if(n_temps==10){
    n_subjects<-1
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
  stopifnot("number of observations must be 20 for the experimental design with 1 TMT-10plex" = n_obs==20)

  #modified_matrix[,3]<-seq(1:n_temps)

  biovar<-(re_var)/(1/icc-1)

  #select  the abundances at the specified temperatures and replicate the fitted values for each TMT 10-plex and each condition
  if(n_temps==1&!design=="OnePot"){#in this case the input has 1 row and 40 cols, each column being a separate subject
    #take the temperature index from the input
    template<-template|>dplyr::group_by(Condition,TechRepMixture)|>dplyr::group_split()
    template<-lapply(template,function(x) x|>dplyr::mutate(Abundance=rep(x$Abundance[t_range],n_replicates/n_conditions),
                                                           temperature=rep(x$temperature[t_range],n_replicates/n_conditions)))|>
      dplyr::bind_rows()


    # create an empty matrix to store modified responses
    modified_matrix <- matrix(nrow = 1, ncol = 20)

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
    result$Condition<-result$treatment
    result$Channel<-unique(template$Channel)[t_range]
    result$TechRepMixture<-1
  }else if(n_temps==2){#If this is a 2 temperature design
    n_temps<-2
    n_subjects<-10
    # create an empty matrix to store modified responses
    modified_matrix <- matrix(NA, nrow = 20, ncol = n_conditions + 2)
    template_temp<-template%>%dplyr::arrange(temperature)
    temperatures <- unique(template_temp$temperature)[t_range]#get temperatures from the data
    for (i in 1:n_conditions) {
      for (j in 1:n_subjects) {
        sample_u <- rnorm(1, mean = 0, sd = sqrt(biovar))#simulate subject effect
        for (k in 1:n_temps) {
          current_responses <- template$Abundance[template$Condition == i][t_range[k]]
          sample_e <- rnorm(1, mean = 0, sd = sqrt(re_var))#simulate random error
          Y_true <- as.numeric(current_responses)
          Y_observed <- Y_true + sample_u + sample_e#add biological variance and error term
          #for a separate function, unlog Y_observed, average Y_obs for all temps, log
          modified_matrix[(j - 1)*n_temps + k , i] <- Y_observed
          modified_matrix[(j - 1)*n_temps + k , n_conditions + 1] <- temperatures[k]#set temperature
          modified_matrix[(j - 1)*n_temps + k , n_conditions + 2] <- j#set Subject
        }
      }
    }
    colnames(modified_matrix) <- c(seq(n_conditions), "temperature","Subject")
    modified_df <- as.data.frame(modified_matrix)

    result<-modified_df|>as.data.frame()|>tidyr::pivot_longer(cols=colnames(modified_df)[colnames(modified_df)!="Subject"&colnames(modified_df)!="temperature"],
                                                              values_to="Abundance",
                                                              names_to="Condition")|>
      dplyr::group_by(Subject,Condition)|>
      dplyr::mutate(Condition=ifelse(Condition==1,"vehicle","treated"))|>
      dplyr::group_split()
    if(design=="OnePot"){
      YonePot_cr<-purrr::map2(result,seq(length(result)),function(x,y) {
        X_obs_crt<-x|>dplyr::mutate(Abundance=2^Abundance)#unlog
        XonePotcr<-X_obs_crt|>dplyr::mutate(Abundance=mean(Abundance,na.rm=T))#average
        YonePot_cr<-XonePotcr|>dplyr::mutate(Abundance=log2(Abundance))#log back
        #collapse temperature after pooling
        YonePot_cr<-YonePot_cr|>
          dplyr::select(-temperature)|>
          dplyr::distinct()|>
          dplyr::mutate(BioReplicate=y,
                        Subject=y,
                        Mixture=Condition)

      })
      result<-dplyr::bind_rows(YonePot_cr)
      #return(YonePot_cr)
    }else if(design=="TPP"){
      result<-dplyr::bind_rows(result)
      #generate temperature values selected
      result<-result|>dplyr::mutate(BioReplicate=Subject,
                                    Mixture=Condition,
                                    Run=Condition)|>
        dplyr::arrange(Condition)
      result$TechRepMixture<-1
      Channel_subject_mapping<-data.frame(Channel=rep(set_temps(10,c(1,2,3,4,5,6,7,8,9,10))$Channel,2),
                                          Subject=unique(result$Subject))|>dplyr::distinct()
      result<-result|>dplyr::inner_join(Channel_subject_mapping,relationship="many-to-many")|>
        dplyr::mutate(Mixture=ifelse(Condition==1,"vehicle","treated"))
      #return(result)
    }

    #generate temperature values selected
    result<-dplyr::bind_rows(result)
    result$Run<-"OnePot"
    result$TechRepMixture<-1

    result$Mixture<-result$Condition
    result$temperature<-mean(unique(modified_df$temperature))
    result<-result|>dplyr::arrange(Condition)
    result$Channel<-as.character(rep(set_temps(10,seq(1,10))$Channel,2))
    result<-dplyr::bind_rows(result)
    return(result)
  }else if(n_temps>1){# if the design has more than one temperature, add variance components

    result<-add_variance_components(template,biovar,re_var,n_conditions,n_subjects,n_temps,t_range)

  }

  return(result)
}
