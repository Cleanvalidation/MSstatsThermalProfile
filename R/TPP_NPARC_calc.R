#Run TPP or NPARC calculation pipeline, renames data from MSstats to TPP format and then uses method = TPP or NPARC to apply the pipeline
TPP_NPARC_calc<-function(protein,method="NPARC",DF=5,temps=NA,CARRIER=TRUE,returnModels=FALSE,filters=FALSE,NORM=TRUE,pipeline=FALSE){
  #rename data to MSStatsTMT format
  start=proc.time()
  protein_All_Data<-list(ProteinLevelData=suppressWarnings(MSstats_format_TO_TPP(protein,temps=temps,CARRIER=CARRIER)))
  end=proc.time()
  print(paste0("Renamed data to match",as.character(method)," format in  ",as.numeric(signif((end-start)[1],2))," seconds"))
  if(method=="NPARC"){
    filters=FALSE#NPARC does not use R^2 or plateaufor a spline model
    stopifnot(any(DF %in% c(3,4,5,6,7)))
    #simulate data for spline DF=5
    start=proc.time()
    protein_NPARC<-runTPP_splineFtest(protein_All_Data$ProteinLevelData,DF=DF,models=returnModels)
    end=proc.time()

    print(paste0("TPP spline F-test results calculated in ",as.numeric(signif((end-start)[1],2)),"seconds"))
    protein_NPARC$unmoderatedFp_val<-tryCatch(1-pf(protein_NPARC$F_statistic, protein_NPARC$df1,protein_NPARC$df2),error=function(cond){return(NA)})
    return(protein_NPARC)

  }else if(method=="TPP"){
    stopifnot(any(DF %in% c(3,4,5,6,7)))
    #simulate data for sigmoid DF=5
    start=proc.time()
    protein_TPP<-runTPP_sigmoid(protein_All_Data$ProteinLevelData,DF=DF,filters=filters,NORM=NORM,pipeline=pipeline)
    end=proc.time()
    print(paste0("TPP results calculated in ",as.numeric(signif((end-start)[1],2))," seconds"))
    return(protein_TPP)
  }
}

