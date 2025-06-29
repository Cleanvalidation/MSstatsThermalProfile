add_variation_shift_2plex <- function(template, icc, n_conditions, n_replicates, re_var, t_range, design = "TPP") {
  library(dplyr)
  library(tidyr)
  library(stringr)

  n_temps <- length(t_range)


  template$treatment <- tolower(template$treatment)
  template <- template %>% mutate(Condition = ifelse(treatment == "vehicle", "1", "2"))


  if (design == "TPP") {
    # Each temperature needs a channel. 10 temps * 2 subjects = 20 channels max.
    n_subjects<-switch(as.character(n_temps),
                       "1" = {
                         20
                       },
                       "2" = {
                         10
                       },
                       "5" = {
                         4
                       },
                       "10" = {
                         2
                       },
                       {
                         stop("Unsupported number of temperatures. Must be 1, 2, 5, or 10.")
                       })
  } else if (design == "OnePot" && n_temps == 2) {

    n_subjects <- (n_conditions * n_replicates)
    n_obs <- n_subjects

  } else if (design == "OnePot" && n_temps == 5)  {
    n_subjects <- (n_conditions * n_replicates)
    n_obs <- n_subjects
  }




  stopifnot("Number of observations must be 20 for the experimental design with 2 TMT-10plexes" = n_obs == 20)

  biovar <- re_var / (1 / icc - 1)

  if (n_temps == 1 && design != "OnePot") {
    template <- template %>% group_by(Condition, TechRepMixture) %>% group_split()
    template <- lapply(template, function(x) x %>%
                         mutate(Abundance = rep(x$Abundance[t_range], n_replicates / n_conditions),
                                temperature = rep(x$temperature[t_range], n_replicates / n_conditions))) %>%
      bind_rows()

    modified_matrix <- matrix(nrow = 1, ncol = 20)

    for (i in 1:n_conditions) {
      for (j in 1:n_subjects) {
        sample_u <- rnorm(1, 0, sqrt(biovar))
        current_responses <- template$Abundance[template$Subject == unique(template$Subject)[j]]
        sample_e <- rnorm(1, 0, sqrt(re_var))
        Y_true <- as.numeric(current_responses[1])
        Y_observed <- Y_true + sample_u + sample_e
        modified_matrix[1, j + (i - 1) * n_subjects] <- Y_observed
      }
    }

    colnames(modified_matrix) <- paste0(rep(t_range, n_subjects * n_conditions), "_", seq_len(n_subjects * n_conditions))
    modified_df <- as.data.frame(modified_matrix)

    result <- pivot_longer(modified_df, cols = everything(), names_to = "Subject", values_to = "Abundance") %>%
      mutate(
        Condition = ifelse(str_detect(Subject, "vehicle"), 1, 2),
        TechRepMixture = 1,
        Temperature = rep(t_range, nrow(.)),
        BioReplicate = Subject,
        Mixture = Subject,
        Run = Subject,
        treatment = ifelse(Condition == "1", "vehicle", "treatment"),
        temperature = str_remove(str_extract(Subject, "[[:digit:]]+_"), "_"),
        Condition = treatment,
        Channel = rep(unique(template$Channel)[seq_along(t_range)], each = nrow(.) / length(t_range))
      )
    return(result)

  } else if (n_temps == 2) {
    n_subjects <- n_replicates
    modified_matrix <- matrix(NA, nrow = 20,ncol=n_conditions + 2)
    template_temp <- template %>% arrange(temperature)
    temperatures <- sort(unique(template_temp$temperature))[t_range]

    for (i in 1:n_conditions) {
      for (j in 1:n_subjects) {
        sample_u <- rnorm(1, 0, sqrt(biovar))

        for (k in 1:n_temps) {
          current_responses <- template$Abundance[template$Condition == i][t_range[k]]
          sample_e <- rnorm(1, 0, sqrt(re_var))
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

    result <- bind_rows(result)
    result$Run <- "OnePot"
    result$TechRepMixture <- 1
    result$Mixture <- result$Condition
    result$temperature <- mean(unique(modified_df$temperature))
    result <- result %>% arrange(Condition)
    result$Channel <- as.character(rep(set_temps(10, 1:10)$Channel, length.out = nrow(result)))
    return(result)
  }

  if (n_temps > 1) {
    result <- add_variance_components(template, biovar, re_var, n_conditions, n_subjects, n_temps, t_range)
  }

  return(result)
}
