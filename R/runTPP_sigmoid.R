#Run TPP sigmoids on normalized data
runTPP_sigmoid<-function(Data,DF=future::availableCores()*0.5,NORM=FALSE,filters=FALSE,pipeline=FALSE){

  if(any(names(Data$TPPdata[[1]])=="uniqueID")){
    Data$TPPdata <- furrr::future_map(Data$TPPdata,function(x) x|>
                                        dplyr::select(-dplyr::starts_with("uniqueID"),
                                                      -dplyr::starts_with("variable"),
                                                      -dplyr::starts_with("value"))|>
                                        as.data.frame())
  }
  Data$TPPdata<-Data$TPPdata|>purrr::map(function(x) x|>dplyr::filter(!stringr::str_detect(gene_name,";")))
  #Data$TPPdata<-lapply(Data$TPPdata,function(x) x[duplicated(x$gene_name),])
  Data$TPPdata<-lapply(Data$TPPdata,function(x) tibble::column_to_rownames(x,"gene_name"))
  Data$TPPdata<-lapply(Data$TPPdata,function(x) x|>dplyr::mutate(gene_name=rownames(x)))
  tpptrData<-TPP::tpptrImport(configTable = Data$TPPconfig,data=Data$TPPdata)
  TRreqs<-TPP::tpptrDefaultNormReqs()
  if(!isTRUE(filters)){


    TRreqs$fcRequirements$thresholdLower<-rep(-1,3)
    TRreqs$fcRequirements$thresholdUpper<- rep(Inf,3)
  }
  if(isFALSE(pipeline)){
  if(isTRUE(NORM)){

    normData<-TPP::tpptrNormalize(data=tpptrData,normReqs = TRreqs)
    Fit_sig <- TPP::tpptrCurveFit(data=normData$normData,
                                  resultPath=getwd(),
                                  nCores = DF,
                                  verbose=FALSE,
                                  maxAttempts = 20,
                                  doPlot=FALSE)
  }else{
    Fit_sig <- TPP::tpptrCurveFit(data=tpptrData,
                                  resultPath=getwd(),
                                  nCores = DF,
                                  verbose=FALSE,
                                  maxAttempts = 20,
                                  doPlot=TRUE)
  }
  if(isTRUE(filters)){
    pval_results<-TPP::tpptrAnalyzeMeltingCurves(Fit_sig,pValFilter = list(minR2 = 0.8, maxPlateau = 0.3))

  }else{
    pval_results<-TPP::tpptrAnalyzeMeltingCurves(Fit_sig,pValFilter = list(minR2 = 0, maxPlateau = 1))

  }
  }else{#if running the entire pipeline
   pval_results<-analyzeTPPTR(
     Data$TPPconfig,
     data = Data$TPPdata,
     resultPath = NULL,
     methods = "meltcurvefit",
     idVar = "gene_name",
     fcStr = "rel_fc_",
     ciStr = NULL,
     naStrs = c("NA", "n/d", "NaN", "<NA>"),
     qualColName = "qupm",
     normalize = FALSE,
     normReqs = TRreqs,
     ggplotTheme = tppDefaultTheme(),
     nCores = "max",
     startPars = c(Pl = 0, a = 550, b = 10),
     splineDF = 5,
     maxAttempts = 500,
     plotCurves = TRUE,
     fixedReference = NULL,
     pValMethod = "robustZ",
     pValFilter = list(minR2 = 0.8, maxPlateau = 0.3),
     pValParams = list(binWidth = 300),
     verbose = FALSE,
     xlsxExport = TRUE
   )

  }
  return(pval_results)
}
