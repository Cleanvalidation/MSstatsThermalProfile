add_variation_shift_2plex <- function(template, icc, n_conditions, n_replicates, re_var, t_range, design = "TPP") {
  library(dplyr)
  library(tidyr)
  library(stringr)

  n_temps <- length(t_range)


  template$treatment <- tolower(template$treatment)
  template <- template %>% mutate(Condition = ifelse(treatment == "vehicle", "1", "2"))


  if (design == "TPP") {
    # Each temperature needs a channel. 10 temps * 2 subjects = 20 channels max.
    if (n_temps %in% c(1, 2, 5, 10)) {
      n_subjects <- 20 / (n_temps * n_conditions)
      if (n_subjects != floor(n_subjects)) {
        stop("Number of subjects is not an integer. Reduce n_temps or check design.")
      }
    } else {
      stop("Please define n_temps as 1, 2, 5, or 10")
    }
  } else if (design == "OnePot" && n_temps == 2) {

    n_subjects <- 20 / (n_conditions * n_replicates)
    n_obs <- n_subjects

  } else if (design == "OnePot" && n_temps == 5)  {
    n_subjects <- 20 / (n_conditions * n_replicates)
    n_obs <- n_subjects * n_conditions * n_replicates  # or just n_obs <- 20
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
    modified_matrix <- matrix(NA, nrow = 20, ncol = n_conditions + 2)
    template_temp <- template %>% arrange(temperature)
    temperatures <- unique(template_temp$temperature)[t_range]

    for (i in 1:n_conditions) {
      for (j in 1:n_subjects) {
        sample_u <- rnorm(1, 0, sqrt(biovar))
        for (k in 1:n_temps) {
          current_responses <- template$Abundance[template$Condition == i][t_range[k]]
          sample_e <- rnorm(1, 0, sqrt(re_var))
          Y_true <- as.numeric(current_responses)
          Y_observed <- Y_true + sample_u + sample_e
          modified_matrix[(j - 1) * n_temps + k, i] <- Y_observed
          modified_matrix[(j - 1) * n_temps + k, n_conditions + 1] <- temperatures[k]
          modified_matrix[(j - 1) * n_temps + k, n_conditions + 2] <- j
        }
      }
    }

    colnames(modified_matrix) <- c(seq_len(n_conditions), "temperature", "Subject")
    modified_df <- as.data.frame(modified_matrix)

    result <- modified_df %>%
      pivot_longer(cols = -c(Subject, temperature), values_to = "Abundance", names_to = "Condition") %>%
      group_by(Subject, Condition) %>%
      mutate(Condition = ifelse(Condition == 1, "vehicle", "treated")) %>%
      group_split()

    if (design == "OnePot") {
      result <- purrr::map2(result, seq_along(result), function(x, y) {
        x %>% mutate(Abundance = 2^Abundance) %>%
          mutate(Abundance = log2(mean(Abundance, na.rm = TRUE))) %>%
          select(-temperature) %>%
          distinct() %>%
          mutate(BioReplicate = y, Subject = y, Mixture = Condition)
      }) %>% bind_rows()
    } else if (design == "TPP") {
      result <- bind_rows(result) %>%
        mutate(BioReplicate = Subject, Mixture = Condition, Run = Condition, TechRepMixture = 1)
      Channel_subject_mapping <- data.frame(
        Subject = unique(result$Subject),
        Channel = rep(set_temps(10, 1:10)$Channel, length.out = length(unique(result$Subject)))
      )
      result <- result %>% inner_join(Channel_subject_mapping, by = "Subject") %>%
        mutate(Mixture = ifelse(Condition == 1, "vehicle", "treated"))
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
