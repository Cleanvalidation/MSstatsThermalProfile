test_that("TPP example runs", {

  data(hdacTR_smallExample, package = "TPP")
  tpptrInput<-TPP::tpptrImport(configTable = hdacTR_config,data=hdacTR_data)

  TRreqs<-TPP::tpptrDefaultNormReqs()
  TRreqs$fcRequirements<-TRreqs$fcRequirements[1,]
  TRreqs$fcRequirements$thresholdLower<-c(-1)
  TRreqs$fcRequirements$thresholdUpper<-c(Inf)
  tpptrData<-TPP::tpptrNormalize(data=tpptrInput,normReqs = TRreqs)
  testthat::expect_s4_class(tpptrData$normData$Vehicle_1,"ExpressionSet")
  normTPP<-tpptrData$normData
  norm_MSStats_format<-TPPnorm_rename(input=tpptrInput,df.temps=Channel2Temps_HumanData)
  testthat::expect_equal(c("Experiment","Protein","Accession","TechRepMixture","Condition","Abundance","Channel","Mixture","BioReplicate","Run","temperature"),names(norm_MSStats_format))

})
