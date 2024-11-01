sigmoid_errors=function(proteinDF,proteinName,cond="vehicle",n_protein){

  all=proteinDF|>dplyr::group_by(Protein)
  # Define the sigmoid model
  Cetsa = function(p,k,m,t) (1-p)/(1 + exp(-k*(1/t - 1/m))) + p

  proteinCond<-all[stringr::str_detect(all$Condition,cond),]|>dplyr::mutate(Protein=as.character(Protein))


  proteinCond=proteinCond|>
    dplyr::group_split()|>
    purrr::map(function(x) x|>dplyr::mutate(fit=list(tryCatch(stats::nls(Abundance ~ Cetsa(p, k, m, temperature),
                                                         data = x,
                                                         start = c(p = 0.05, k = 985, m = 50),
                                                         control = nlme::nlmeControl(msMaxIter = 5000,tol = 1e-5),
                                                         na.action = na.omit),error=function(cond){return(NA)}))))


  #Keep the number of proteins that converged from the sigmoid model

  proteinCond_keep<-proteinCond|>purrr::keep(function(x)any(!is.na(x$fit[[1]])))

  #print the % of proteins that fit the sigmoid
  print(paste0(length(unique(dplyr::bind_rows(proteinCond_keep)$Protein))," out of ",length(unique(dplyr::bind_rows(proteinCond)$Protein))," proteins fit the sigmoid model"))
  #filter out missing values from each protein profile
  res<-purrr::map(proteinCond_keep,function(x) x|>
                         dplyr::select(uniqueID,temperature,fit,Abundance,Subject)|>
                    dplyr::filter(!is.na(Abundance)))
  #generate residuals for each protein profile
  temp_res<-purrr::map(proteinCond_keep,function(x) data.frame(
    residuals=as.vector(
      t(resid(x$fit[[1]]))))|>
    dplyr::mutate(Subject=x$Subject[!is.na(x$Abundance)],
                  uniqueID=x$uniqueID[1]))
  LMM<-purrr::map(temp_res,function(x){LMM=print(
      lme4::VarCorr(lme4::lmer(residuals ~ 1 + (1 | Subject),data=x),comp="Variance"))|>
      as.data.frame()
    df<-data.frame(r_Var=LMM$vcov[2]^2,
                   s_Var=LMM$vcov[1]^2,
                   icc=LMM$vcov[1]^2/(LMM$vcov[2]^2+LMM$vcov[1]^2),
                   uniqueID=x$uniqueID[1])
    return(df)
  })
  #join original data to the variance components and ICC results
  data<-purrr::map2(
  proteinCond_keep,LMM,function(x,y)x|>dplyr::inner_join(y))

  #Make a histogram of the ICC for the actual dataset used in the Sims
  ICC<-lapply(data,function(x) x[1,])|>dplyr::bind_rows()
  hist(ICC$icc,breaks=10,xlab = "ICC",col="#D95F0E",cex.lab = 1.5,cex.axis=1.5,main="TPP proc. ICC")

  #calculate the median residual error from LMM
  r_error<-sqrt(median(ICC$r_Var))

  #define a s_Var range from ICC values 0.05, 0.1,0.2, 0.3, 0.4 ICC= s_var/(r_var+s_var)
  range<-c(0.05,0.1,0.2,0.4,0.6)
  s_var<-range*(r_error)^2/(1-range)
  re_var<-(r_error)^2
  #select a template protein
  template<-data|>purrr::keep(function(x) any(x$uniqueID==proteinName))

  #make profiles with partition residual error and subject error for each protein template
  Sims<-lapply(template,function(x)
    make_profiles(x,n_protein,s_Var=s_var,r_Var=re_var))
  Sims<-dplyr::bind_rows(Sims)

  #Benchmarks on TPP
  resultTPP<-TPP_NPARC_calc(Sims,method="NPARC",temps=temps,CARRIER=FALSE,returnModels = FALSE,filters=FALSE,NORM=FALSE)
  hist(resultTPP$p_NPARC_unmod,breaks=20,xlab = "p-values",col="#D95F0E",cex.lab = 1.5,cex.axis=1.5,main="TPP F-test unmod")
  #Benchmarks on NPARC
  resultNPARC<-NPARC::NPARCfit(x = Sims$temperature,
                               y = Sims$Abundance,
                               id = Sims$Protein,
                               groupsNull = NULL,
                               groupsAlt = as.character(Sims$Condition),
                               returnModels = FALSE)
  testNPARC<-NPARC::NPARCtest(resultNPARC$metrics,dfType="theoretical")
  hist(testNPARC$pVal,breaks=20,xlab = "p-values",col="#FEC44F",cex.lab = 1.5,cex.axis=1.5,main="NPARC F-test")

  #Benchrmaks on MSstats

  Sims$Run<-Sims$Spectrum.File<-paste0(Sims$Condition,"_",Sims$TechRepMixture)
  Sims$Condition<-paste0(Sims$Channel,"_",Sims$Condition)
  Sims$Condition<-ifelse(Sims$temperature==min(Sims$temperature,na.rm=T),"Norm",Sims$Condition)
  Sims$Mixture<-Sims$Subject
  Sims$BioReplicate<-Sims$Condition
  Sims$temperature<-as.character(Sims$temperature)
  dataMSstat<-list(ProteinLevelData=Sims)
  comparison<-make_contrast_matrix(dataMSstat,temps=c("53","56","59"))

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
  hist(ATE_MSstats$ComparisonResult$pvalue,breaks=20,xlab = "p-values",col="#2C7FB8",cex.lab = 1.5,cex.axis=1.5,main="MSstatsTMT")

   #join icc results
  testNPARC<-testNPARC|>dplyr::mutate(uniqueID=id,Protein=id)|>dplyr::inner_join(Sims)|>dplyr::mutate(icc=formatC(icc, digits = 2))
  resultTPP<-resultTPP|>dplyr::mutate(Protein=uniqueID)|>dplyr::inner_join(Sims|>dplyr::select(uniqueID,Protein,Subject,Condition,TechRepMixture,bio_var,error_var,icc)|>dplyr::distinct())|>dplyr::mutate(icc=formatC(icc,digits = 2))
  resultMSstats<-ATE_MSstats$ComparisonResult|>dplyr::mutate(uniqueID=Protein)|>
    dplyr::inner_join(Sims|>dplyr::select(uniqueID,Protein,Subject,Condition,TechRepMixture,bio_var,error_var,icc)|>dplyr::distinct())|>
                        dplyr::select(Protein,uniqueID,icc,pvalue)|>
                        dplyr::distinct()|>dplyr::mutate(icc=formatC(icc,digits = 2))


  testNPARC<-testNPARC|>dplyr::select(uniqueID,pVal,icc,fStat)|>dplyr::distinct()
  resultTPP<-resultTPP|>dplyr::select(uniqueID,p_NPARC,p_NPARC_unmod,icc,F_statistic)|>dplyr::distinct()
  ggplot(testNPARC,mapping=aes(x=pVal))+
    geom_histogram(color="black",fill="#FEC44F",bins=20)+facet_wrap(~icc,nrow=1)+xlab("unmoderated p-value")+theme(text=element_text(size=20))

  ggplot(resultTPP,mapping=aes(x=p_NPARC_unmod))+
    geom_histogram(color="black",fill="#D95F0E",bins=20)+facet_wrap(~icc,nrow=1)+xlab("unmoderated p-value")+theme(text=element_text(size=20))

  ggplot(resultTPP,mapping=aes(x=p_NPARC))+
    geom_histogram(color="black",fill="#D95F0E",bins=20)+facet_wrap(~icc,nrow=1)+xlab("moderated p-value")+theme(text=element_text(size=20))

  ggplot(resultMSstats,mapping=aes(x=pvalue))+
    geom_histogram(color="black",fill="#2C7FB8",bins=20)+facet_wrap(~icc,nrow=1)+xlab("unmoderated p-value")+theme(text=element_text(size=20))

  FPNPARC<-testNPARC$uniqueID[testNPARC$pVal<0.05]
  FPTPP<-resultTPP$uniqueID[resultTPP$p_NPARC<0.05]
#
#   All<-Sims|>
#     dplyr::mutate(bin_bvar=cut(bio_var,4),
#                   id=Protein)|>
#     dplyr::inner_join(testNPARC)|>
#
#     dplyr::filter(id%in% FPNPARC)|>
#     dplyr::group_by(bin_bvar)|>
#     dplyr::group_split()|>
#     purrr::map(function(x) x|>dplyr::filter(id %in% unique(x$id)[1:5]))|>
#     dplyr::bind_rows()
#   #plot false positive profiles
#   FPNPARC_plots<-ggplot(All,
#                         mapping=aes(x=temperature,y=Abundance,color=Condition,shape=Subject))+
#     geom_point(size=2.5)+
#     facet_wrap(~Protein)

  return(list(TPP=resultTPP,NPARC=testNPARC))
}

