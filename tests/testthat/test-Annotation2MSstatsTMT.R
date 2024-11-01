test_that("Annotation file function is ready for MSStatsTMT", {
  data("PSMsample_HumanData", package = "MSstatsThermalProfiler")
  pd_annotation<-Annotation2MSstatsTMT(PSMsample_HumanData,solvent="DMSO",temps=Channel2Temps_HumanData,reference="126")

  processed.input <- MSstatsTMT::PDtoMSstatsTMTFormat(PSMsample_HumanData,
                                                      pd_annotation,
                                                      which.proteinid = "Master Protein Accessions",
                                                      useNumProteinsColumn = TRUE,
                                                      useUniquePeptide = FALSE,
                                                      rmPSM_withfewMea_withinRun = TRUE,
                                                      rmProtein_with1Feature = TRUE,
                                                      summaryforMultipleRows = max,
                                                      use_log_file = FALSE,
                                                      append=FALSE,
                                                      verbose=TRUE)
  testthat::expect_equal(names(processed.input),c( "ProteinName","PeptideSequence","Charge","PSM","Mixture","TechRepMixture","Run","Channel","BioReplicate","Condition","Intensity"))
})
