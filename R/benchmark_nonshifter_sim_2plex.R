benchmark_nonshifter_sim_2plex <- function(msstats_icc_output, templateProtein, n_sims, t_range = seq(7, 7), design = "hybrid") {
  # Load required packages
  if (!requireNamespace("MSstatsTMT", quietly = TRUE)) stop("The MSstatsTMT package is required but not installed.")
  if (!requireNamespace("MSstats", quietly = TRUE)) stop("The MSstats package is required but not installed.")
  if (!requireNamespace("MSstatsConvert", quietly = TRUE)) stop("The MSstatsConvert package is required but not installed.")

  # Get protein data with variance
  all_proteins <- msstats_icc_output$df_with_variance
  template_MSstats <- list(ProteinLevelData = dplyr::filter(all_proteins, Protein %in% templateProtein))

  # Rename for MSstats
  template_MSstats$ProteinLevelData$LogIntensities <- template_MSstats$ProteinLevelData$Abundance
  template_MSstats$ProteinLevelData$GROUP <- template_MSstats$ProteinLevelData$Condition
  template_MSstats$ProteinLevelData$SUBJECT <- template_MSstats$ProteinLevelData$BioReplicate

  # MSstats quantification
  template_MSstats <- MSstats::quantification(template_MSstats, type = "Group", use_log_file = FALSE)
  template_MSstats <- tidyr::pivot_longer(template_MSstats,
                                          cols = setdiff(names(template_MSstats), "Protein"),
                                          names_to = "Condition",
                                          values_to = "Abundance")

  # Get template simulation data
  template_simulation <- template_MSstats|>
    dplyr::arrange(temperature)

  # Match channel-temperature mapping
  temps <- all_proteins |>
    dplyr::filter(Protein %in% templateProtein) |>
    dplyr::select(Channel, temperature) |>
    dplyr::arrange(temperature)|>
    dplyr::distinct()

  template_simulation <- template_simulation |>
    tidyr::separate(Condition, into = c("Channel", "treatment")) |>
    dplyr::mutate(Condition = paste0(Channel, "_", treatment),
                  BioReplicate = Condition) |>
    dplyr::inner_join(temps, by = "Channel") |>
    dplyr::group_by(temperature) |>
    dplyr::mutate(Abundance = mean(Abundance, na.rm = TRUE)) |>
    dplyr::ungroup()

  # Plot template
  png(filename = "template_MsstatsTMTproc_nonshifter.png",
      width = 12, height = 12, units = "in", pointsize = 12,
      res = 600, type = "cairo")
  Template <- ggplot2::ggplot(template_simulation, ggplot2::aes(x = treatment, y = Abundance, color = treatment)) +
    ggplot2::geom_point(size = 4) +
    #ggplot2::geom_step(linewidth = 1) +
    #ggplot2::scale_x_continuous("Temperature", breaks = unique(template_simulation$temperature)) +
    ggplot2::labs(title = "Simulation template: non shifter", x = "Temperature", y = "Log of protein abundances") +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 24),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  print(Template)
  dev.off()
  saveRDS(Template, "templateMSstatsTMTproc_nonshifter.RDS")

  # Setup simulation
  icc_range <- c(0.05, 0.2, 0.4, 0.6, 0.8)
  re_var <- median(all_proteins$sigma_re^2, na.rm = TRUE)

  # Define conditions and temperatures
  n_conditions <- 2
  n_temps <- length(t_range)

  # Calculate valid number of replicates
  max_channels <- 20
  n_replicates <- floor(max_channels / (n_conditions * n_temps))
  n_obs <- n_replicates * n_conditions * n_temps

  if (n_obs != max_channels) {
    stop(paste("Number of observations must be", max_channels,
               "for the experimental design with 2 TMT-10plexes. Got", n_obs, "instead."))
  }

  template_simulation$Subject <- interaction(template_simulation$treatment, template_simulation$temperature)

  # Format template for simulation
  df <- template_simulation |>
    dplyr::select(Abundance, temperature, treatment, Channel) |>
    dplyr::arrange(temperature)

  result <- list()

  for (i in seq_len(n_sims)) {
    for (j in icc_range) {
      Sims <- add_variation_shift_2plex(df, j, n_conditions, n_replicates, re_var, t_range, design = design)
      Sims <- dplyr::bind_rows(Sims)
      Sims$Protein <- paste0("Sim_", i, "_icc_", j)
      result <- append(result, list(Sims))
    }
  }

  result <- dplyr::bind_rows(result)
  return(result)
}
