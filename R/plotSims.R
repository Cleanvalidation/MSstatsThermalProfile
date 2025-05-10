#' @importFrom stringr str_extract


plotSims<-function(df,labels=NA,processing="",temps=c("53.8","57.1","60.4")){
  if (!requireNamespace("MSstatsTMT", quietly = TRUE)) {
    stop("The MSstatsTMT package is required but not installed.")
  }
  if (!requireNamespace("MSstats", quietly = TRUE)) {
    stop("The MSstats package is required but not installed.")
  }
  if (!requireNamespace("MSstatsConvert", quietly = TRUE)) {
    stop("The MSstatsConvert package is required but not installed.")
  }

  if(!length(unique(df$Condition))==2){
    df$Condition<-stringr::str_extract(stringr::str_extract_all(stringr::str_to_lower(df$Condition),"_[[:lower:]]+"),"[[:lower:]]+")
    df$Condition<-ifelse(df$Condition=="character","Norm",df$Condition)
    df$Condition<-df$treatment<-ifelse(df$Condition=="vehicle","vehicle","treated")
  }else{
    df$Condition<-df$treatment<-ifelse(df$treatment=="vehicle","vehicle","treated")

  }

  dataMSstat<-list(ProteinLevelData=df)

  comparison<-make_contrast_matrix_all(dataMSstat,temps=NA)

  ATE_MSstats<-MSstatsTMT::groupComparisonTMT(
    dataMSstat,
    contrast.matrix = comparison,
    moderated = FALSE,
    adj.method = "BH",
    remove_norm_channel = TRUE,
    remove_empty_channel = TRUE,
    save_fitted_models = TRUE,
    use_log_file = FALSE,
    append = FALSE,
    verbose = TRUE,
    log_file_path = NULL
  )
  ATE_MSstats$ComparisonResult$ICC<-stringr::str_extract(ATE_MSstats$ComparisonResult$Protein,"icc_[:digit:].[[:digit:]]+")
  ATE_MSstats$ComparisonResult$ICC<-paste0("% of bio var = ",100*as.numeric(stringr::str_extract(ATE_MSstats$ComparisonResult$ICC,"[:digit:].[[:digit:]]+")))
  ATE_MSstats$ComparisonResult<-ATE_MSstats$ComparisonResult[stringr::str_detect(ATE_MSstats$ComparisonResult$ICC,"5|40"),]
  ATE_MSstats$ComparisonResult<-ATE_MSstats$ComparisonResult|>dplyr::group_by(ICC)|>dplyr::mutate(Sens=100*sum(pvalue<0.001)/length(unique(Protein)))
  df<-df|>
    as.data.frame()|>
    dplyr::group_by(Protein,Run,Mixture)|>
    dplyr::group_split()|>
    lapply(function(x) x|>
             dplyr::mutate(Abundance=2^Abundance))


  df<-lapply(df,function(y) y|>as.data.frame()|>dplyr::mutate(Abundance=Abundance/Abundance[10]))|>dplyr::bind_rows()

  ATE_MSstats$ComparisonResult<-ATE_MSstats$ComparisonResult[stringr::str_detect(ATE_MSstats$ComparisonResult$ICC,"5|40|80"),]
  Msstat<-ggplot(ATE_MSstats$ComparisonResult,mapping=aes(x=pvalue))+
    geom_histogram(fill="#2C7FB8",color="black")+facet_wrap(~ICC,nrow=1)+
    theme(text=element_text(size=15))+
    scale_x_continuous(n.breaks=3)+ylim(0,250)+xlab("unmoderated p-value")
  #Benchmarks on TPP
  proteins1<-result|>dplyr::mutate(Condition=stringr::str_extract(Condition,"[[:lower:]]+"),
                                   shape=ifelse(temperature==min(temperature,na.rm=TRUE),"reference","included"),
  )|>dplyr::distinct()
  proteins1$Subject<-ifelse(stringr::str_ends(proteins1$Subject,"_1"),
                            paste0(proteins1$Condition,"_",1),
                            paste0(proteins1$Condition,"_",2))
  #unlog abundances
  proteins1<-proteins1|>
    as.data.frame()|>
    dplyr::group_by(Protein,Run,Mixture)|>
    dplyr::group_split()|>
    lapply(function(x) x|>
             dplyr::mutate(Abundance=2^Abundance))
  if(length(colnames(comparison)[comparison!=0])>=10){
    proteins1<-dplyr::bind_rows(proteins1)|>dplyr::group_by(Protein)|>dplyr::group_split()
  }
  #scale abundances
  proteins1<-lapply(proteins1,function(y) y|>
                      as.data.frame()|>
                      dplyr::mutate(Abundance=
                                      Abundance/Abundance[y$shape=="reference"]))|>
    dplyr::bind_rows()
  TPPdata<-proteins1|>
    dplyr::mutate(uniqueID=Protein,
                  Condition=treatment,
                  Subject=as.character(Subject))
  resultTPP<-TPP_NPARC_calc(TPPdata,method="NPARC",temps=temps,CARRIER=FALSE,returnModels = FALSE,filters=FALSE,NORM=FALSE)
  resultTPP<-resultTPP[stringr::str_detect(resultTPP$ICC,"5|40|80"),]
  Tpp<-ggplot(resultTPP,mapping=aes(x=p_NPARC))+
    geom_histogram(fill="#D95F0E",color="black")+facet_wrap(~ICC,nrow=1)+
    theme(text=element_text(size=15))+xlab("moderated p-value")+ylim(0,250)
  Tpp_unmod<-ggplot(resultTPP,mapping=aes(x=p_NPARC_unmod))+
    geom_histogram(fill="#D95F0E",color="black")+facet_wrap(~ICC,nrow=1)+
    theme(text=element_text(size=15))+xlab("unmoderated p-value")+ylim(0,250)
  Tpp_manual_F<-ggplot(resultTPP,mapping=aes(x=unmoderatedFp_val))+
    geom_histogram(fill="#D95F0E",color="black")+facet_wrap(~ICC,nrow=1)+
    theme(text=element_text(size=15))+xlab("unmoderated p-value")+ylim(0,250)
  #Benchmarks on NPARC
  resultNPARC<-NPARC::NPARCfit(x = TPPdata$temperature,
                               y = TPPdata$Abundance,
                               id = TPPdata$Protein,
                               groupsNull = NULL,
                               groupsAlt = as.character(TPPdata$treatment),
                               returnModels = FALSE)
  testNPARC<-NPARC::NPARCtest(resultNPARC$metrics,dfType="theoretical")
  testNPARCemp<-NPARC::NPARCtest(resultNPARC$metrics,dfType="empirical")
  testNPARCemp<-testNPARCemp[stringr::str_detect(testNPARCemp$ICC,"5|40|80"),]
  Nparc<-ggplot(testNPARCemp,mapping=aes(x=pVal))+
    geom_histogram(fill="#FEC44F",color="black")+facet_wrap(~ICC,nrow=1)+
    theme(text=element_text(size=15))+xlab("scaled df p-value")+ylim(0,250)
  #Benchmarks on SCAM
  ATE_SCAM<-compute_pvalues_ATE_RE_F(TPPdata)
  Scam<-ATE_SCAM
  Scam<-Scam[stringr::str_detect(Scam$ICC,"5|40|80"),]
  Scam<-ggplot(Scam,mapping=aes(x=p.value))+
    geom_histogram(fill="#07fff8",color="black")+facet_wrap(~ICC,nrow=1)+
    theme(text=element_text(size=15))+xlab("unmoderated p-value")+ylim(0,250)
  #plot profiles
  resultTPP<-resultTPP|>
    dplyr::filter(uniqueID %in% One_prot_ICC$Protein)|>
    dplyr::select(uniqueID,p_NPARC_unmod)|>
    dplyr::distinct()|>
    dplyr::mutate(Protein=as.character(uniqueID),
                  p_TPP_unmod=p_NPARC_unmod)|>
    dplyr::select(-p_NPARC_unmod)
  resultNPARC<-testNPARC|>
    dplyr::filter(id %in% One_prot_ICC$Protein)|>
    dplyr::select(id,pVal)|>
    dplyr::distinct()|>
    dplyr::mutate(Protein=as.character(id),
                  p_NPARC_unmod=pVal)|>
    dplyr::select(-pVal)
  resultSCAM<-ATE_SCAM|>
    dplyr::filter(Accession %in% One_prot_ICC$Protein)|>
    dplyr::select(Accession,p.value)|>
    dplyr::distinct()|>
    dplyr::mutate(Protein=as.character(Accession),
                  p_SCAM_unmod=p.value)|>
    dplyr::select(-p.value)
  resultMSstatsTMT<-ATE_MSstats$ComparisonResult|>
    dplyr::select(Protein,pvalue)|>
    dplyr::distinct()|>
    dplyr::mutate(Protein=as.character(Protein),
                  p_MSstat_unmod=pvalue)|>
    dplyr::select(-pvalue)
  ICC_list<-One_prot_ICC|>
    dplyr::group_by(Protein)|>
    dplyr::inner_join(resultTPP)|>
    dplyr::inner_join(resultNPARC)|>
    dplyr::inner_join(resultSCAM)|>
    dplyr::inner_join(resultMSstatsTMT)|>dplyr::group_split()
  ICC_vocab<-list(ICC=c("A","B","C","D"),
                  ICC=c("E","F","G","H"),
                  ICC=c("I","J","K","L"),
                  ICC=c("M","N","O","P"),
                  ICC=c("Q","R","S","T"))  #these are the legend panels
  plotList_ICC<-purrr::map2(ICC_list,ICC_vocab,function(x,y)
    plotSims(x,labels=y,processing="MSstatsTMT",temps=c("53.8","57.1","60.4")))
  return(P_fits)

}

