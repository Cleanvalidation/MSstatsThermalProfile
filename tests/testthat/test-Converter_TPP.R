test_that("Converter works", {
  data("PSMsample_HumanData", package = "MSstatsThermalProfiler")
  names(PSMsample_HumanData)<-stringr::str_replace_all(names(PSMsample_HumanData),"[[:punct:]]+| ",".")
  #data<-PSMsample_HumanData[,!stringr::str_detect(names(PSMsample_HumanData),"Abundance..131C|Abundance.131C")]
  Channel2Temps_HumanData$Channel<-stringr::str_replace(Channel2Temps_HumanData$Channel,"131N","131")
  names(PSMsample_HumanData)<-stringr::str_replace(names(PSMsample_HumanData),"131N","131")
  pd_annotation<-Annotation2MSstatsTMT(PSMsample_HumanData,solvent="DMSO",temps=Channel2Temps_HumanData,reference="126",CARRIER=TRUE)


  #PSMs to peptide groups
  processed.input <- MSstatsTMT::PDtoMSstatsTMTFormat(pd_annotation$df,
                                                      pd_annotation$annotation,
                                                      which.proteinid = "Master Protein Accessions",
                                                      useNumProteinsColumn = FALSE,
                                                      useUniquePeptide = FALSE,
                                                      rmPSM_withfewMea_withinRun = TRUE,
                                                      rmProtein_with1Feature = TRUE,
                                                      summaryforMultipleRows = max,
                                                      use_log_file = FALSE,
                                                      log_file_path=NULL,
                                                      append=FALSE,
                                                      verbose=TRUE)

  #peptide groups to Proteins
  summarised.proteins <- MSstatsTMT::proteinSummarization(data = processed.input,
                                                          method = "msstats",
                                                          global_norm = FALSE,
                                                          reference_norm = TRUE,
                                                          MBimpute=FALSE,
                                                          use_log_file=FALSE,
                                                          append=FALSE,
                                                          verbose=FALSE,
                                                          remove_norm_channel=FALSE)
  temps<-set_temps(10,c(37.3, 40.6, 43.9, 47.2, 50.5, 53.8, 57.1, 60.4, 64, 67))
  x<-Converter_TPP(summarised.proteins$ProteinLevelData,CARRIER=TRUE)
  expect_length(x,16)
})
