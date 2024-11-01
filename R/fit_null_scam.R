fit_null_scam = function(df){
  model = scam::scam(Abundance ~ s(temperature,bs="mpd",k=5),data=df)
  return(model)
}
