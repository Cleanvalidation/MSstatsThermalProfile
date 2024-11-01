find_deltaTm = function(simdata_shifted){
  params = simdata_shifted$params

  n = nrow(params)

  params$Tm_tmt = NaN
  params$Tm_ctrl = NaN

  for(i in 1:n){
    p_Tmt = params$p_Tmt[i]
    k_Tmt = params$k_Tmt[i]
    m_Tmt = params$m_Tmt[i]

    p_Ctrl = params$p_Ctrl[i]
    k_Ctrl = params$k_Ctrl[i]
    m_Ctrl = params$m_Ctrl[i]

    Tm_tmt = uniroot(function(t) Cetsa(p_Tmt,k_Tmt,m_Tmt,t) - 0.5,interval=c(30,70))$root
    Tm_ctrl = uniroot(function(t) Cetsa(p_Ctrl,k_Ctrl,m_Ctrl,t) - 0.5,interval=c(30,70))$root

    params$Tm_tmt[i] = Tm_tmt
    params$Tm_ctrl[i] = Tm_ctrl
  }

  params = params |> dplyr::mutate(delta_Tm = Tm_tmt - Tm_ctrl)
  return(params)
}
