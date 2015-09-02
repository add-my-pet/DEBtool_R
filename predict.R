#-------------------
#------------------- Translated from matlab code of Kooijman, see http://www.bio.vu.nl/thb/deb/deblab/add_my_pet_old/html/Oncorhynchus_mykiss.html
#-------------------



predict_Oncorhynchus_mykiss = function (par, data, tW) {
  

  # par: (21)-vector with parameters (see below)
  # data:(r_d,1)-matrix with zero-variate data (not some elements are used)
  # tW: matrix with uni-variate data with time since birth, wet weight
  # Edata:(r_d,1)-matrix with expected values for data(:,2)
  # EW: vectors with expected values for tW(:,2)
      
  ## unpack par
  T_ref  = par[1] # K, temp for which rate pars are given 
  T_A    = par[2] # K, Arrhenius temp
  f      = par[3] # -, scaled functional response
  z      = par[4] # -, zoom factor
  del_M  = par[5] # -, shape coefficient to convert vol-length to physical length
  F_m    = par[6] # l/d.cm^2, {F_m}, max spec searching rate
  kap_X  = par[7] # -, digestion efficiency of food to reserve
  v      = par[8] # cm/d, energy conductance
  kap    = par[9] # -, allocation fraction to soma
  kap_R  = par[10]# -, reproduction efficiency
  p_M    = par[11]# J/d.cm^3, [p_M], vol-specific somatic maintenance
  p_T    = par[12]# J/d.cm^2, {p_T}, surface-specific som maintenance
  k_J    = par[13]# 1/d, maturity maint rate coefficient
  E_G    = par[14]# J/cm^3, [E_G], spec cost for structure
  E_Hh   = par[15]# J, E_H^j, maturity at hatch
  E_Hb   = par[16]# J, E_H^b, maturity at birth
  E_Hj   = par[17]# J, E_H^j, maturity at metamorphosis
  E_Hp   = par[18]# J, E_H^p, maturity at puberty
  h_a    = par[19]# 1/d^2, Weibull aging acceleration
  s_G    = par[20]# -, Gompertz stress coefficient
  f_tW   = par[21]# -, scaled functional response for tL-data
  W_0    = par[22]# g, wet weight at 0
  
  d_V = dwm[2,1] 
  d_E = dwm[3,1]
  w_V = dwm[2,2] 
  w_E = dwm[3,2]
  mu_V = dwm[2,3] 
  mu_E = dwm[3,3]
  M_V = d_V/ w_V                  # mol/cm^3, [M_V], volume-specific mass of structure
  kap_G = M_V * mu_V/ E_G         # -, growth efficiency
  
  # Selected copy-paste from parscomp & statistics
  p_Am = z * p_M/ kap             # J/d.cm^2, {p_Am} spec assimilation flux
  k_M = p_M/ E_G                  # 1/d, somatic maintenance rate coefficient
  k = k_J/ k_M                    # -, maintenance ratio
  # p_Xm = p_Am/ kap_X            # J/d.cm^2, max spec feeding power
  
  y_V_E = mu_E * M_V/ E_G         # mol/mol, yield of structure on reserve
  y_E_V = 1/ y_V_E                # mol/mol, yield of reserve on structure
  
  E_m = p_Am/ v                   # J/cm^3, reserve capacity [E_m]
  m_Em = y_E_V * E_m/ E_G         # mol/mol, reserve capacity 
  g = E_G/ (kap* E_m)             # -, energy investment ratio
  w = m_Em * w_E/ w_V             # -, contribution of reserve to weight
  
  L_m = v/ k_M/ g                 # cm, maximum length
  L_T = p_T/ p_M                  # cm, heating length (also applies to osmotic work)
  l_T = L_T/ L_m                  # -, scaled heating length
  
  J_E_Am = p_Am/ mu_E             # mol/d.cm^2, {J_EAm}, max surface-spec assimilation flux
  # maturity at hatch
  M_Hh = E_Hh/ mu_E               # mol, maturity at hatch
  U_Hh = M_Hh/ J_E_Am             # cm^2 d, scaled maturity at hatch
  V_Hh = U_Hh/ (1 - kap)          # cm^2 d, scaled maturity at hatch
  v_Hh = V_Hh * g^2 * k_M^3/ v^2  # -, scaled maturity at hatch
  u_Hh = U_Hh * g^2 * k_M^3/ v^2  # -, scaled maturity at hatch
  # maturity at birth
  M_Hb = E_Hb/ mu_E               # mol, maturity at birth
  U_Hb = M_Hb/ J_E_Am             # cm^2 d, scaled maturity at birth
  V_Hb = U_Hb/ (1 - kap)          # cm^2 d, scaled maturity at birth
  v_Hb = V_Hb * g^2 * k_M^3/ v^2  # -, scaled maturity at birth
  u_Hb = U_Hb * g^2 * k_M^3/ v^2  # -, scaled maturity at birth
  # maturity at metamorphosis
  M_Hj = E_Hj/ mu_E  		   # mol, maturity at metamorposis
  U_Hj = M_Hj/ J_E_Am             # cm^2 d, scaled maturity at metamorphosis 
  V_Hj = U_Hj/ (1 - kap)          # cm^2 d, scaled maturity at metamorposis
  v_Hj = V_Hj * g^2 * k_M^3/ v^2  # -, scaled maturity at metamorphosis
  u_Hj = U_Hj * g^2 * k_M^3/ v^2  # -, scaled maturity at metamorphosis 
  # maturity at puberty
  M_Hp = E_Hp/ mu_E               # mol, maturity at puberty
  U_Hp = M_Hp/ J_E_Am             # cm^2 d, scaled maturity at puberty 
  V_Hp = U_Hp/ (1 - kap)          # cm^2 d, scaled maturity at puberty
  v_Hp = V_Hp * g^2 * k_M^3/ v^2  # -, scaled maturity at puberty
  u_Hp = U_Hp * g^2 * k_M^3/ v^2  # -, scaled maturity at puberty  
  
  ## zero-variate data
  
  # temperature correct for all data
  # hatch
  TC = tempcorr(data[1,1], T_ref, T_A)    # -, Temperature Correction factor
  pars_lb = c(g,k,v_Hh)        # compose parameter vector
  tblbinfo = get_tb(pars_lb, f)    # -, scaled age and length at birth at f
  t_h=tblbinfo[1]
  l_h=tblbinfo[2]
  info=tblbinfo[3]
  # if info ~= 1 # numerical procedure failed
  # fprintf('warning: invalid parameter value combination for get_tb \n')
  # end
  L_h = L_m * l_h     # cm, structural length at hatch at f
  a_h = t_h/ k_M 
  aT_h = a_h/ TC     # d, age at hatch 
  
  # birth
  TC = tempcorr(data[2,1], T_ref, T_A)    # -, Temperature Correction factor
  pars_lj =c(g,k,l_T,v_Hb,v_Hj,v_Hp)
  
  tjtptbljlplblirjrBinfo = get_tj(pars_lj, f)
  t_j=tjtptbljlplblirjrBinfo[1]
  t_p=tjtptbljlplblirjrBinfo[2]
  t_b=tjtptbljlplblirjrBinfo[3]
  l_j=tjtptbljlplblirjrBinfo[4]
  l_p=tjtptbljlplblirjrBinfo[5]
  l_b=tjtptbljlplblirjrBinfo[6]
  l_i=tjtptbljlplblirjrBinfo[7]
  rho_j=tjtptbljlplblirjrBinfo[8]
  rho_b=tjtptbljlplblirjrBinfo[9]
  info=tjtptbljlplblirjrBinfo[10]
  
  # if info ~= 1 # numerical procedure failed
  # fprintf('warning: invalid parameter value combination for get_tb \n')
  # end
  L_b = L_m * l_b                         # cm, structural length at birth at f
  Lw_b = L_b/ del_M                       # cm, physical length at birth at f
  a_b = t_b/ k_M                      # d, age at birth at f
  aT_b = a_b/ TC
  
  # metamorphosis
  L_j = l_j * L_m 
  Lw_j = L_j/ del_M
  
  # puberty 
  TC = tempcorr(data[3,1], T_ref, T_A)    # -, Temperature Correction factor
  L_p = l_p * L_m 
  Lw_p = L_p/ del_M      # cm, structural, physical length at puberty
  aT_p = t_p/ k_M/ TC                     # d, age at puberty
  
  # ultimate size
  L_i = L_m * l_i 
  Lw_i = L_i/ del_M      # cm, ultimate structural, physical length
  Ww_i = L_i^3 * (1 + f * w)              # g, ultimate wet weight
  
  # reproduction
  TC = tempcorr(data[9,1], T_ref, T_A)    # -, Temperature Correction factor
  pars_R = c(kap,kap_R,g,k_J,k_M,L_T,v,U_Hb,U_Hp,L_b,L_j,L_p) # compose parameter vector
  RT_i = TC * reprod_rate_metam(L_i, f, pars_R)[1] # ultimate reproduction rate
  
  # life span
  TC = tempcorr(data[4,1], T_ref, T_A)    # -, Temperature Correction factor
  pars_tm = c(g,l_T, h_a/k_M^2, s_G)     # compose parameter vector
  t_m = get_tm_s(pars_tm, f, l_b, l_p)    # -, scaled mean life span
  aT_m = t_m/ k_M/ TC                     # d, mean life span
  
  # pack output for zero-variate data
  Edata = c(aT_h, aT_b, aT_p,  aT_m, Lw_b, Lw_p, Lw_i, Ww_i, RT_i, # real data
           v, kap, kap_R, p_M, p_T, k_J, kap_G)                 # pseudo data
  
  
  ## uni-variate data
  # t-Ww-data
  TC = tempcorr(273 + 8.5, T_ref, T_A)             # -, Temperature Correction factor
  
  tjtptbljlplblirjrBinfo = get_tj(pars_lj, f_tW)
  t_j=tjtptbljlplblirjrBinfo[1]
  t_p=tjtptbljlplblirjrBinfo[2]
  t_b=tjtptbljlplblirjrBinfo[3]
  l_j=tjtptbljlplblirjrBinfo[4]
  l_p=tjtptbljlplblirjrBinfo[5]
  l_b=tjtptbljlplblirjrBinfo[6]
  l_i=tjtptbljlplblirjrBinfo[7]
  rho_j=tjtptbljlplblirjrBinfo[8]
  rho_b=tjtptbljlplblirjrBinfo[9]
  info=tjtptbljlplblirjrBinfo[10]
  
  rT_B = TC * rho_b * k_M # 1/d, von Bert, exponential growth rate
  rT_j = TC * rho_j * k_M 
  aT_b = t_b/ k_M/ TC 
  aT_j = t_j/ k_M/ TC
  L_b = l_b * L_m
  L_j = l_j * L_m 
  L_i = l_i * L_m
  L_0 = (W_0/ (1 + f_tW * w))^(1/3) 
  aT_0 = aT_b + log(L_0/ L_b) * 3/ rT_j
  t_j = aT_j - aT_0 # time since zero at metamorphosis
  t_bj = tW[which(tW[,1] < t_j),1] # select times between birth & metamorphosis
  L_bj = L_0 * exp(t_bj * rT_j/3) # exponential growth as V1-morph
  t_ji = tW[which(tW[,1] >= t_j),1] # selects times after metamorphosis
  L_ji = L_i - (L_i - L_j) * exp( - rT_B * (t_ji - t_j)) # cm, expected length at time
  L = c(L_bj, L_ji) # catenate lengths
  EW = L^3 * (1 + f_tW * w)
  
  return(list(Edata, EW))
}

