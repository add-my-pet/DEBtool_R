## get_tm_s
# Obtains scaled mean age at death fro short growth periods

##
#  function [tm, Sb, Sp, info] = 
get_tm_s = function(p, f, lb, lp) {
  # created 2009/02/21 by Bas Kooijman, modified 2014/03/17, 2015/01/18
  
  ## Syntax
  # [tm, Sb, Sp, info] = <../get_tm_s.m *get_tm_s*>(p, f, lb, lp)
  
  ## Description
  # Obtains scaled mean age at death assuming a short growth period relative to the life span
  # Divide the result by the somatic maintenance rate coefficient to arrive at the mean age at death. 
  # The variant get_tm_foetus does the same in case of foetal development.
  # If the input parameter vector has only 4 elements (for [g, lT, ha/ kM2, sG]), 
  #   it skips the calulation of the survival probability at birth and puberty.   
  #
  # Input
  #
  # * p: 4 or 7-vector with parameters: [g lT ha sG] or [g k lT vHb vHp ha SG]
  # * f: optional scalar with scaled reserve density at birth (default eb = 1)
  # * lb: optional scalar with scaled length at birth (default: lb is obtained from get_lb)
  # * lp: optional scalar with scaled length at puberty
  #  
  # Output
  #
  # * tm: scalar with scaled mean life span
  # * Sb: scalar with survival probability at birth (if length p = 7)
  # * Sp: scalar with survival prabability at puberty (if length p = 7)
  # * info: indicator equals 1 if successful, 0 otherwise
  
  ## Remarks
  # Theory is given in comments on DEB3 Section 6.1.1. 
  # See <get_tm.html *get_tm*> for the general case of long growth period relative to life span
  
  ## Example of use
  # get_tm_s([.5, .1, .1, .01, .2, .1, .01])
  
  if (!exists('f')) {
    f = 1
  } else if (is.na(f)) {
    f = 1;
  }
  
  if (length(p) >= 7){
    #  unpack pars
    g   = p[1] # energy investment ratio
    #k   = p[2] # k_J/ k_M, ratio of maturity and somatic maintenance rate coeff
    lT  = p[3] # scaled heating length {p_T}/[p_M]Lm
    #vHb = p[4] # v_H^b = U_H^b g^2 kM^3/ (1 - kap) v^2; U_B^b = M_H^b/ {J_EAm}
    #vHp = p[5] # v_H^p = U_H^p g^2 kM^3/ (1 - kap) v^2; U_B^p = M_H^p/ {J_EAm}
    ha  = p[6] # h_a/ k_M^2, scaled Weibull aging acceleration
    sG  = p[7] # Gompertz stress coefficient
  } else if (length(p) == 4){
    #  unpack pars
    g   = p[1] # energy investment ratio
    lT  = p[2] # scaled heating length {p_T}/[p_M]Lm
    ha  = p[3] # h_a/ k_M^2, scaled Weibull aging acceleration
    sG  = p[4] # Gompertz stress coefficient
  }
  
  if (abs(sG) < 1e-10) {
    sG = 1e-10
  }
  
  li = f - lT
  hW3 = ha * f * g/ 6/ li
  hW = hW3^(1/3) # scaled Weibull aging rate
  hG = sG * f * g * li^2
  hG3 = hG^3     # scaled Gompertz aging rate
  tG = hG/ hW
  tG3 = hG3/ hW3             # scaled Gompertz aging rate
  
  info_lp = 1
  if (length(p) >= 7) {
    if (!exists('lp') & !exists('lb')) {
      lplbinfo = get_lp(p, f)
      lp=lplbinfo[1]
      lb=lplbinfo[2]
      info_lp=lplbinfo[3]
    } else if (!exist('lp') || is.na('lp')){
      lplbinfo = get_lp(p, f, lb)
      lp=lplbinfo[1]
      lb=lplbinfo[2]
      info_lp=lplbinfo[3]
    }
  
  
  # get scaled age at birth, puberty: tb, tp
    tblbinfo_tb = get_tb(p, f, lb)
    tb=tblbinfo_tb[1]
    lb=tblbinfo_tb[2]
    info_tb=tblbinfo_tb[3]
    
    irB = 3 * (1 + f/ g)
    tp = tb + irB * log((li - lb)/ (li - lp))
    hGtb = hG * tb
    Sb = exp((1 - exp(hGtb) + hGtb + hGtb^2/ 2) * 6/ tG3)
    hGtp = hG * tp
    Sp = exp((1 - exp(hGtp) + hGtp + hGtp^2/ 2) * 6/ tG3) 
    if (info_lp == 1 && info_tb == 1){
      info = 1
    } else{
      info = 0
    }
  } else {  # length(p) == 4
    Sb = NaN
    Sp = NaN
    info = 1
  }
  
  if (abs(sG) < 1e-10) {
    tm = gamma(4/3)/ hW
    tm_tail = 0
  } else if (hG > 0) {
    tm = 10/ hG # upper boundary for main integration of S(t)
    library(pracma)
    tm_tail = expint(exp(tm * hG) * 6/ tG3)/ hG
    tm = quad(fnget_tm_s, 0, tm * hW, tol = 1.0e-12, trace = FALSE,tG)/ hW + tm_tail
  } else {# hG < 0
    tm = -1e4/ hG # upper boundary for main integration of S(t)
    hw = hW * sqrt( - 3/ tG) # scaled hW
    tm_tail = sqrt(pi)* erfc(tm * hw)/ 2/ hw
    tm = quad(fnget_tm_s, 0, tm * hW, tol = 1.0e-12, trace = FALSE, tG)/ hW + tm_tail
  }
}
  
  # subfunction
  
  #function S = 
fnget_tm_s=function(t, tG){
  # modified 2010/02/25
  # called by get_tm_s for life span at short growth periods
  # integrate ageing surv prob over scaled age
  # t: age * hW 
  # S: ageing survival prob
  
  hGt = tG * t # age * hG
  if (tG > 0) {
    a= matrix(hGt) %*% t(matrix((1/(4:500))))
    S = exp(-(1 + apply((t(apply(a, 1, cumprod))),1, sum)) * t^3)             #apply((t(apply(a, 1, cumprod))),2, sum)?
  } else {# tG < 0
    S = exp((1 + hGt + hGt^2/ 2  - exp(hGt)) * 6/ tG^3)
  }
  return(S) 
}