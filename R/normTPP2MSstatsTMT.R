#normalizes and renames output from MSstats_format_TO_TPP
normTPP2MSstatsTMT<-function(df,df.temps=Channel2Temps_HumanData){
  tpptrInput<-TPP::tpptrImport(configTable = df$TPPconfig,data=df$TPPdata)

  TRreqs<-TPP::tpptrDefaultNormReqs()
  TRreqs$fcRequirements<-TRreqs$fcRequirements[1,]
  TRreqs$fcRequirements$thresholdLower<-c(-1)
  TRreqs$fcRequirements$thresholdUpper<-c(Inf)
  tpptrData<-TPP::tpptrNormalize(data=tpptrInput,normReqs = TRreqs)

  normTPP<-tpptrData$normData
  norm_MSStats_format<-TPPnorm_rename(normTPP,df.temps)
  return(norm_MSStats_format)
}
