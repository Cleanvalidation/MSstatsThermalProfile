compute_spline_ATE_singleprotein = function(data){
  model = lm(Abundance ~ splines::ns(temperature,df=5) + Condition,data=data)
  result = broom::tidy(model) |> dplyr::filter(term=="Conditiontreated") |> dplyr::rename(pvalue=p.value,SE=std.error)
  return(result)
}

compute_spline_ATE = function(data){
  gdf = data |> dplyr::group_by(Protein)
  keys = (gdf |> dplyr::group_keys())$Protein
  result = gdf |> dplyr::group_split() |> purrr::map_dfr(.f=function(x) compute_spline_ATE_singleprotein(x))
  result$Protein = keys
  return(result)
}

compute_spline_Ftest_singleprotein = function(data){
  null = lm(Abundance ~ splines::ns(temperature,df=5),data=data)
  alt = lm(Abundance ~ splines::ns(temperature,df=5)*Condition,data=data)
  result = broom::tidy(model) |> dplyr::filter(term=="Conditiontreated") |> dplyr::rename(pvalue=p.value,SE=std.error)
  return(result)
}
