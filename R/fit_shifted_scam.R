fit_shifted_scam = function(df){
  model = scam::scam(Abundance ~ s(temperature,bs="mpd",by=Condition) + Condition,data=df)
  return(model)
}
