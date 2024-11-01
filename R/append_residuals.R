#append residuals of a single sigmoid fit across conditions to the protein data structure
append_residuals<-function(all_proteins){
  #Get all protein names
  proteinNames<-unique(all_proteins$Protein)
  #count n proteins
  n_proteins<-length(proteinNames)
  #add residual column to all_proteins
  all_proteins$residuals<-rep(NA,nrow(all_proteins))
  #define model
  Cetsa = function(p,k,m,t) (1-p)/(1 + exp(-k*(1/t - 1/m))) + p
  #select protein i
  for (i in 1:n_proteins){
    #print(i)
    protein_i<-all_proteins[all_proteins$Protein == proteinNames[i],]
    #print(protein_i)
    #fit nls
    model_fit<-NA
    model_fit<-tryCatch(stats::nls(Abundance ~ Cetsa(p, k, m, temperature),
                        data = protein_i,
                        start = c(p = 0.05, k = 985, m = 50),
                        control = nlme::nlmeControl(msMaxIter = 5000,tol = 1e-5),
                        na.action = na.omit),error=function(x){
                          return(NA)
                        })
    #print(model_fit)
    if(!is.na(model_fit)){
      all_proteins$residuals[all_proteins$Protein == proteinNames[i]]<-residuals(model_fit)
    }
    return(all_proteins)
  }

}
