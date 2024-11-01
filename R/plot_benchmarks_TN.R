plot_benchmarks_TN<-function(data,resultTPP,resultNPARC,resultSCAM,resultMSstatsTMT){
  resultMSstatsTMT$ComparisonResult<-resultMSstatsTMT$ComparisonResult|>dplyr::mutate(ICC=NA)
  resultMSstatsTMT$ComparisonResult$ICC<-stringr::str_extract(resultMSstatsTMT$ComparisonResult$Protein,"icc_[:digit:].[[:digit:]]+")
  resultMSstatsTMT$ComparisonResult$ICC<-paste0("% of bio var = ",100*as.numeric(stringr::str_extract(resultMSstatsTMT$ComparisonResult$ICC,"[:digit:].[[:digit:]]+")))
  resultMSstatsTMT$ComparisonResult<-resultMSstatsTMT$ComparisonResult[stringr::str_detect(resultMSstatsTMT$ComparisonResult$ICC,"5|40|80"),]
  resultMSstatsTMT$ComparisonResult$ICC<-factor(resultMSstatsTMT$ComparisonResult$ICC,levels=c("% of bio var = 5","% of bio var = 40","% of bio var = 80"))

  Msstat<-ggplot(resultMSstatsTMT$ComparisonResult,mapping=aes(x=pvalue))+
    geom_histogram(fill="#2C7FB8",color="black")+facet_wrap(~ICC,nrow=3)+
    theme(text=element_text(size=15))+
    scale_x_continuous(n.breaks=3)+ylim(0,150)+xlab("unmoderated p-value")+
    ggtitle("MSstatsTMT proc. MSstatsTMT stat model.")+
    scale_x_continuous(n.breaks=8)+ylab("protein count")

  #Benchmarks on TPP
  resultTPP$ICC<-stringr::str_extract(resultTPP$uniqueID,"icc_[:digit:].[[:digit:]]+")
  resultTPP$ICC<-paste0("% of bio var = ",100*as.numeric(stringr::str_extract(resultTPP$ICC,"[:digit:].[[:digit:]]+")))
  resultTPP<-resultTPP[stringr::str_detect(resultTPP$ICC,"5|40|80"),]
  resultTPP$ICC<-factor(resultTPP$ICC,levels=c("% of bio var = 5","% of bio var = 40","% of bio var = 80"))
  Tpp<-ggplot(resultTPP,mapping=aes(x=p_NPARC))+
    geom_histogram(fill="#D95F0E",color="black")+facet_wrap(~ICC,nrow=3)+
    theme(text=element_text(size=15))+xlab("moderated p-value")+ylim(0,150)+ggtitle("MSstatsTMT proc. TPP stat model.")+ scale_x_continuous(n.breaks=8)+ylab("protein count")

  Tpp_unmod<-ggplot(resultTPP,mapping=aes(x=p_NPARC_unmod))+
    geom_histogram(fill="#D95F0E",color="black")+facet_wrap(~ICC,nrow=3)+
    theme(text=element_text(size=15))+xlab("unmoderated p-value")+ylim(0,150)

  Tpp_manual_F<-ggplot(resultTPP,mapping=aes(x=unmoderatedFp_val))+
    geom_histogram(fill="#D95F0E",color="black")+facet_wrap(~ICC,nrow=3)+
    theme(text=element_text(size=15))+xlab("unmoderated p-value")+ylim(0,150)
  #Benchmarks on NPARC
  resultNPARC$ICC<-stringr::str_extract(resultNPARC$id,"icc_[:digit:].[[:digit:]]+")
  resultNPARC$ICC<-paste0("% of bio var = ",100*as.numeric(stringr::str_extract(resultNPARC$ICC,"[:digit:].[[:digit:]]+")))
  resultNPARC<-resultNPARC[stringr::str_detect(resultNPARC$ICC,"5|40|80"),]
  resultNPARC$ICC<-factor(resultNPARC$ICC,levels=c("% of bio var = 5","% of bio var = 40","% of bio var = 80"))
  Nparc<-ggplot(resultNPARC,mapping=aes(x=pVal))+
    geom_histogram(fill="#FEC44F",color="black")+facet_wrap(~ICC,nrow=3)+
    theme(text=element_text(size=15))+xlab("scaled df p-value")+ylim(0,150)+ggtitle("MSstatsTMT proc. NPARC stat model.")+ scale_x_continuous(n.breaks=8)+ylab("protein count")

  #Benchmarks on SCAM
  resultSCAM$ICC<-stringr::str_extract(resultSCAM$Accession,"icc_[:digit:].[[:digit:]]+")
  resultSCAM$ICC<-paste0("% of bio var = ",100*as.numeric(stringr::str_extract(resultSCAM$ICC,"[:digit:].[[:digit:]]+")))
  resultSCAM<-resultSCAM[stringr::str_detect(resultSCAM$ICC,"5|40|80"),]

  resultSCAM<-resultSCAM[stringr::str_detect(resultSCAM$ICC,"5|40|80"),]
  resultSCAM$ICC<-factor(resultSCAM$ICC,levels=c("% of bio var = 5","% of bio var = 40","% of bio var = 80"))
  Scam<-ggplot(resultSCAM,mapping=aes(x=p.value))+
    geom_histogram(fill="#07fff8",color="black")+facet_wrap(~ICC,nrow=3)+
    theme(text=element_text(size=15))+xlab("unmoderated p-value")+ylim(0,150)+ggtitle("MSstatsTMT proc. SCAM stat model.")+ scale_x_continuous(n.breaks=8)+ylab("protein count")

  histogram_plot<-ggpubr::ggarrange(Tpp,Nparc,Scam,Msstat, ncol=4, nrow=1,
                                  common.legend = TRUE, legend="none")+theme_bw()+theme(panel.border = element_blank())

  #Define a data frame with one protein sim per ICC value
  result<-data|>dplyr::inner_join(resultMSstatsTMT$ComparisonResult|>as.data.frame()|>dplyr::select(Protein,ICC)|>dplyr::distinct())

 return(histogram_plot)
}
