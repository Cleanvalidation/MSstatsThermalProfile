#Run NPARC implementation F test from TPP package
runTPP_splineFtest<-function(normData,DF=5,models=FALSE){
  normData$TPPdata<-purrr::map(normData$TPPdata,function(x) x |> as.data.frame())
  tpptrData<-TPP::tpptrImport(configTable = normData$TPPconfig,data=normData$TPPdata)


  check_longTable <- TPP::tpptrTidyUpESets(tpptrData)
  #TPP fit splines
  #check if condition column is missing after TPPtrImport
  if(all(is.na(check_longTable$condition))){
    check_longTable$condition<-stringr::str_extract(check_longTable$experiment,"[[:lower:]]+")
  }
  SimSplineFits<- TPP::tpptrFitSplines(data = check_longTable,
                                        factorsH1 = c("condition"),
                                        splineDF = DF,
                                        nCores = 6,
                                        returnModels = models)


  #compute F-test
  SimFtest<-TPP::tpptrFTest(SimSplineFits)

  # plotSplines<-TPP::tpptrPlotSplines(
  #   data=check_longTable,
  #   factorsH1 = NULL,
  #   factorsH0 = NULL,
  #   SimSplineFits,
  #   testResults=SimFtest,
  #   resultPath = "~/CS7290/MSstatsThermalProfiler/",
  #   individual = TRUE,
  #   overview = FALSE,
  #   returnPlots = TRUE,
  #   control = list(nCores = "max", maxRank = 500, highlightBelow = 0.05),
  #   maxRank = NULL,
  #   highlightBelow = NULL,
  #   plotIndividual = NULL,
  #   plotAlphabetical = NULL
  # )
  SimFtest$p_NPARC_unmod<-1 - suppressWarnings(pf(SimFtest$F_statistic,
                                                  SimFtest$df1, SimFtest$df2))
  SimFtest$p_adj_NPARC_unmod<-p.adjust(SimFtest$p_NPARC_unmod,
                                       method = "BH")
  SplineFits<-SimSplineFits[SimSplineFits$testHypothesis=="alternative",]

  if(isTRUE(models)){
    SimF<-SimFtest |> dplyr::inner_join(SplineFits)|> dplyr::select(-df1,-df2) %>% dplyr::group_by(uniqueID) %>% dplyr::group_split()
    resid_vs_fit_Splines<-purrr::map2(SimF,SplineFits$fittedModel,function(x,y){ tryCatch({p<-data.frame(Protein=as.factor(x$uniqueID[1]),
                                                                                     pNPARC=x$p_NPARC[1],
                                                                                     padjNPARC=x$p_adj_NPARC[1],
                                                                                     successfulFit=x$successfulFit[1],
                                                                                     sigma=summary(y)$sigma[1])
                                                                                 q<-data.frame(fitted=y$fitted.values,residuals=y$residuals)
                                                                                 r<-cbind(p,q)
                                                                                 return(r)},
                                                                                 error=function(e){
                                                                                   return(NA)})
      })


    resid_vs_fit_Splines<-resid_vs_fit_Splines|>purrr::keep(function(x) any(!is.na(x)))
    resid_vs_fit_Splines<-dplyr::bind_rows(resid_vs_fit_Splines)
    resid_vs_fit_Splines<-resid_vs_fit_Splines|>dplyr::mutate(uniqueID=Protein)|>
      dplyr::inner_join(SimFtest)
    return(resid_vs_fit_Splines)
  }else{

  return(SimFtest)
  }
}
