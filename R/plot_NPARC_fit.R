plot_NPARC_fit = function(data,nparcfit_result,protein,temprange=seq(35,67,1)){
  if(any(names(data)=="id")&!any(names(data)=="Protein")){
    protein_data<-data|>dplyr::mutate(Protein=id)|>dplyr::filter(Protein==protein)
  }else{
    protein_data<-data|>dplyr::filter(Protein==protein)
  }

  protein_params = nparcfit_result$metrics |> dplyr::filter(id==protein[1],modelType=="alternative")

  treated_params = protein_params |> dplyr::filter(!group=="vehicle")
  vehicle_params = protein_params |> dplyr::filter(group=="vehicle")

  cetsa = function(a,b,Pl,x) (1 - Pl)  / (1+exp((b - a/x))) + Pl

  treated_preds = data.frame(Abundance = cetsa(treated_params$a,treated_params$b,
                                               treated_params$pl,temprange))|>
    dplyr::mutate(temperature=temprange,Condition="treated")

  vehicle_preds = data.frame(Abundance = cetsa(vehicle_params$a,vehicle_params$b,
                                               vehicle_params$pl,temprange))|>dplyr::mutate(Condition="vehicle")
  vehicle_preds$temperature<-temprange

  preds_df = dplyr::bind_rows(treated_preds,vehicle_preds)
  P<-ggplot(protein_data,aes(x=temperature,y=Abundance,col=Condition,shape=shape),size=2) + geom_point(size=2.5) +
    geom_line(inherit.aes = FALSE,data=preds_df,aes(x=temperature,y=Abundance,col=Condition),size=1)+ theme(text = element_text(size = 18))
  return(P)
}

#ex
#plot_NPARC_fits(dat,fitresult,"P23458")
