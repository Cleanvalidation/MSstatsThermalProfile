test_that("scam models work", {
  sim=simulate_multiple_shifted_nonsigmoids(runs=10,n=4)
  Sim2MSstat<-Sim2MSstatsTMT(sim)

  if(any(names(Sim2MSstat)=="Protein")){
    Sim2MSstat$Accession=Sim2MSstat$Protein
  }
  if(any(names(Sim2MSstat)=="Condition")&!any(names(Sim2MSstat)=="treatment")){
    Sim2MSstat$treatment=as.factor(Sim2MSstat$Condition)
  }
  if(any(names(Sim2MSstat)=="Abundance")&!any(names(Sim2MSstat)=="I")){
    Sim2MSstat$I=Sim2MSstat$Abundance
  }
  if(any(names(Sim2MSstat)=="Subject")&!any(names(Sim2MSstat)=="sample_id")){
    Sim2MSstat$sample_id=Sim2MSstat$Subject
  }
  if(!any(names(Sim2MSstat)=="sample_name")&any(names(Sim2MSstat)=="Run")){
    Sim2MSstat$sample_name=stringr::str_extract(unique(Sim2MSstat$Run)[length(unique(Sim2MSstat$Run))],"[[:lower:]]+")

  }else{
    Sim2MSstat$sample_name<-"Simulation"
  }
  original_result<-fit_scam_marginal_ATE_RE_F(Sim2MSstat) |> dplyr::bind_rows()
  testthat::expect_equal(names(original_result),c(
                         "Accession","r.sq","dRSS","F_test","F_pvalue","term","contrast","estimate","std.error",
                         "statistic","p.value","conf.low","conf.high"
                         ))
})
