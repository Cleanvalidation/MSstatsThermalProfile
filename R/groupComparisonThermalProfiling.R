#' Finding differentially abundant proteins across conditions in TMT experiment
#'
#' Tests for significant changes in protein abundance across conditions based on a family of linear mixed-effects models in TMT experiment.
#' Experimental design of case-control study (patients are not repeatedly measured) is automatically determined based on proper statistical model.
#'
#' @param data the output of \code{\link{proteinSummarization}} function. It is a list with data frames `FeatureLevelData` and `ProteinLevelData`
#' @param contrast.matrix Comparison between conditions of interests. 1) default is "pairwise", which compare all possible pairs between two conditions. 2) Otherwise, users can specify the comparisons of interest. Based on the levels of conditions, specify 1 or -1 to the conditions of interests and 0 otherwise. The levels of conditions are sorted alphabetically.
#' @param missing_timepoint If "drop", missing time point will be dropped and remaining time-points will be used for comparison.
#' If "replace", missing timepoint will be replaced by available time point provided in the `replacement` parameter.
#' @param replacement Named list of time points that can replace time points provided in the contrast.matrix.
#' @param moderated TRUE will moderate t statistic; FALSE (default) uses ordinary t statistic.
#' @param adj.method adjusted method for multiple comparison. "BH" is default.
#' @param remove_norm_channel TRUE(default) removes "Norm" channels from protein level data.
#' @param remove_empty_channel TRUE(default) removes "Empty" channels from protein level data.
#' @param save_fitted_models logical, if TRUE, fitted models will be added to
#'
#' @return a list that consists of the following elements:
#' (1) ComparisonResult: statistical testing results;
#' (2) FittedModel: the fitted linear models
#'
#' @export
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats anova lm median model.matrix na.omit p.adjust pt
#' @importFrom MSstats MSstatsContrastMatrix
#'
#' @examples
#' data(input.pd)
#' # use protein.summarization() to get protein abundance data
#' quant.pd.msstats = proteinSummarization(input.pd,
#'                                        method="msstats",
#'                                        global_norm=TRUE,
#'                                        reference_norm=TRUE)
#'
#' test.pairwise = groupComparisonTMT(quant.pd.msstats, moderated = TRUE)
#' head(test.pairwise$ComparisonResult)
#'
#' # Only compare condition 0.125 and 1
#' levels(quant.pd.msstats$ProteinLevelData$Condition)
#'
#' # Compare condition 1 and 0.125
#' comparison=matrix(c(-1,0,0,1),nrow=1)
#'
#' # Set the nafmes of each row
#' row.names(comparison)="1-0.125"
#'
#' # Set the column names
#' colnames(comparison)= c("0.125", "0.5", "0.667", "1")
#' test.contrast = groupComparisonTMT(data = quant.pd.msstats,
#' contrast.matrix = comparison,
#' moderated = TRUE)
#' head(test.contrast$ComparisonResult)
#'
groupComparisonThermalProfiling = function(
    data, contrast.matrix = "pairwise",
    missing_timepoint = "drop", replacement = NULL,
    moderated = FALSE, adj.method = "BH",
    remove_norm_channel = TRUE, remove_empty_channel = TRUE,
    save_fitted_models = FALSE,
    use_log_file = TRUE, append = FALSE,
    verbose = TRUE, log_file_path = NULL
){
  MSstatsConvert::MSstatsLogsSettings(
    use_log_file, append, verbose, log_file_path,
    base = "MSstatsTMT_log_groupComparison_", pkg_name = "MSstatsTMT"
  )
  getOption("MSstatsTMTLog")("INFO", "MSstatsTMT - groupComparisonTMT function")
  getOption("MSstatsTMTLog")("INFO", paste("Moderated t-stat :", moderated))
  getOption("MSstatsTMTLog")("INFO", paste("Adjust p-value :", adj.method))
  getOption("MSstatsTMTLog")("INFO", paste("Remove empty channels :",
                                           remove_empty_channel))
  getOption("MSstatsTMTLog")("INFO", paste("Remove normalization channels :",
                                           remove_norm_channel))

  summarized = MSstatsTMT:::MSstatsPrepareForGroupComparisonTMT(data$ProteinLevelData,
                                                                remove_norm_channel,
                                                                remove_empty_channel)
  contrast_matrix = MSstats::MSstatsContrastMatrix(contrast.matrix,
                                                   unique(summarized$Group))
  fitted_models = MSstatsTMT:::MSstatsFitComparisonModelsTMT(summarized)
  FittedModel <- fitted_models$fitted_model
  names(FittedModel) <- fitted_models$protein

  fitted_models = MSstatsTMT:::MSstatsModerateTTest(summarized, fitted_models, moderated)
  testing_results = MSstatsGroupComparisonThermal(fitted_models, contrast_matrix,
                                                  missing_timepoint, replacement)
  testing_results = MSstatsTMT:::MSstatsGroupComparisonOutputTMT(testing_results,
                                                                 adj.method)

  if(save_fitted_models){
    list(ComparisonResult = testing_results,
         ModelQC = NULL,
         FittedModel = FittedModel)
  } else{
    list(ComparisonResult = testing_results,
         ModelQC = NULL,
         FittedModel = NULL)
  }
}


MSstatsGroupComparisonThermal = function(fitted_models, contrast_matrix,
                                         missing_timepoint, replacement) {
  msg = paste0("Testing for ", length(fitted_models) , " proteins:")
  getOption("MSstatsTMTLog")("INFO", msg)
  getOption("MSstatsTMTMsg")("INFO", msg)
  testing_results = vector("list", length(fitted_models))
  pb = txtProgressBar(max = length(fitted_models), style = 3)
  for (i in seq_along(fitted_models)) {
    testing_result = MSstatsTestSingleProteinThermal(fitted_models[[i]],
                                                     contrast_matrix,
                                                     missing_timepoint,
                                                     replacement)
    testing_results[[i]] = testing_result
    setTxtProgressBar(pb, i)
  }
  close(pb)
  testing_results
}


MSstatsTestSingleProteinThermal = function(fitted_model, contrast_matrix,
                                           missing_timepoint, replacement) {
  single_protein = fitted_model[["data"]]
  groups = as.character(unique(single_protein$Group))
  groups = sort(groups)
  coefs = fitted_model[["coefs"]][[1]]
  s2 = fitted_model[["variance"]]
  s2_df = fitted_model[["variance_df"]]
  fit = fitted_model[["fitted_model"]][[1]]
  s2_prior = fitted_model[["variance_prior"]]
  df_prior = fitted_model[["df_prior"]]
  total_df = fitted_model[["total_df"]]
  protein = fitted_model[["protein"]]

  results = vector("list", nrow(contrast_matrix))
  if (!inherits(fit, "try-error")) {

    if (is.infinite(df_prior)) {
      s2_posterior <- s2_prior
    } else{
      s2_posterior = (s2_prior * df_prior + s2 * s2_df) / (df_prior + s2_df)
    }

    for (row_id in seq_len(nrow(contrast_matrix))) {
      contrast = contrast_matrix[row_id, , drop = FALSE]
      positive.groups = colnames(contrast)[contrast > 0]
      negative.groups = colnames(contrast)[contrast < 0]

      if (any(positive.groups %in% groups) &
          any(negative.groups %in% groups)) {

        if (missing_timepoint == "drop") {
          cm = MSstatsTMT:::.makeContrastSingleTMT(fit, contrast, single_protein, coefs)
        } else {
          cm = .makeContrastSingleThermal(fit, contrast, single_protein, coefs,
                                          replacement)
        }
        FC = (cm %*% coefs)[1, 1]
        if (inherits(fit, "lm")) {
          se2.post = diag(t(cm) %*% summary(fit)$cov.unscaled %*% cm) * s2_posterior
          df.post = s2_df + df_prior
        } else {
          # Acknowlege: Tyler Bradshawthis contributed to this part of implementation
          vcov = fit@vcov_beta
          se2 = as.matrix(t(cm) %*% as.matrix(vcov) %*% cm)
          ## calculate posterior variance
          vcov.post = fit@pp$unsc() * s2_posterior
          se2.post = as.matrix(t(cm) %*% as.matrix(vcov.post) %*% cm)
          ## calculate posterior df
          g <- sapply(fit@Jac_list, function(gm) cm %*% gm %*% cm)
          denom <- try(t(g) %*% fit@vcov_varpar %*% g, silent=TRUE)
          if (inherits(denom, "try-error")) {
            df.post = s2_df + df_prior
          } else{
            df.post = 2*(se2)^2/denom + df_prior
          }
        }
        df.post <- pmin(df.post, total_df)

        t = FC / sqrt(se2.post)
        p = 2*pt(-abs(t), df = df.post)
        se = sqrt(se2.post)
        DF = df.post
        if (s2_df == 0) {
          issue = "SingleMeasurePerCondition"
        } else{
          issue = NA
        }
      } else {
        result = MSstatsTMT:::.handleMissingConditionTMT(single_protein, contrast[1, ])
        FC = result[["logFC"]]
        p = NA
        se = NA
        DF = NA
        issue = result[["issue"]]
      }

      results[[row_id]] = list(Protein = protein, Comparison = row.names(contrast),
                               log2FC = FC, pvalue = p, SE = se, DF = DF, issue = issue)
      # This code cannot be extracted to a new function yet due to
      # performance issues with environments
    }
  } else {
    # very few measurements so that the model is unfittable
    for (row_id in 1:nrow(contrast_matrix)) {
      contrast = contrast_matrix[row_id, , drop = FALSE]
      result = MSstatsTMT:::.handleMissingConditionTMT(single_protein, contrast[1, ])
      results[[row_id]] = list(Protein = protein,
                               Comparison = row.names(contrast),
                               log2FC = result$logFC, pvalue = NA, SE = NA,
                               DF = NA, issue = result$issue)
    }
  }
  data.table::rbindlist(results)
}


.makeContrastSingleThermal = function(fit, contrast, single_protein, coefs,
                                      replacement) {
  sub_groups = as.character(unique(single_protein$Group))
  missing_groups = setdiff(colnames(contrast), sub_groups)
  positive = names(contrast)[contrast > 0]
  negative = names(contrast)[contrast < 0]

  replacement_present = lapply(replacement,
                               function(x) {
                                 options = intersect(sub_groups, x)
                                 if (length(options) == 0) {
                                   NULL
                                 } else {
                                   options[1]
                                 }
                               })

  if (!(all(positive %in% sub_groups) & all(negative %in% sub_groups))) {
    for (group in missing_groups) {
      group_value = contrast[, group, drop = T]
      contrast[, group] = 0
      contrast[, replacement_present[[group]]] = group_value
      contrast
    }
  }

  MSstatsTMT:::.makeContrastSingleTMT(fit, contrast, single_protein, coefs)
}

