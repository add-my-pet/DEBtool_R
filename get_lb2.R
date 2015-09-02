## get_lb2
# Obtains scaled length at birth, given the scaled reserve density at birth

##
get_lb2= function(p, eb, lb0=NA){
# created 2007/07/26 by Bas Kooijman; modified 2013/08/19, 2015/01/18

## Syntax
# [lb info] = <../get_lb2.m *get_lb2*> (p, eb, lb0)

## Description
# Obtains scaled length at birth, given the scaled reserve density at birth. 
# Like get_lb, but using the shooting method, rather than Newton Raphson
#
# Input
#
# * p: 3-vector with parameters: g, k, v_H^b (see below)
# * eb: optional scalar with scaled reserve density at birth (default eb = 1)
# * lb0: optional scalar with initial estimate for scaled length at birth (default lb0: lb for k = 1)
#  
# Output
#
# * lb: scalar with scaled length at birth 
# * info: indicator equals 1 if successful, 0 otherwise

## Remarks
# Like <get_lb.html *get_lb*>, but uses a shooting method in 1 variable

#  unpack p
  g = p[1]   # g = [E_G] * v/ kap * {p_Am}, energy investment ratio
  k = p[2]   # k = k_J/ k_M, ratio of maturity and somatic maintenance rate coeff
  vHb = p[3] # v_H^b = U_H^b g^2 kM^3/ (1 - kap) v^2; U_H^b = M_H^b/ {J_EAm}
  
  info = 1
  
  
  if (is.na(lb0)) {
    lb = as.complex(vHb)^(1/ 3) # exact solution for k = 1     
  } else {
    lb = lb0
  }
  if (!exists('eb')){
    eb = 1
  } else if (is.na(eb)){
    eb = 1
  }
  
  
  xb = g/ (eb + g)
  xb3 = xb^(1/3)
  
  xbxb3gvHbk=c(xb, xb3, g, vHb, k)
  lb=Re(lb)
  
  lbflaginfo = fzero(fnget_lb2, lb, maxiter = 100, tol = 10e-8, xbxb3gvHbk)
  #lbflaginfo = uniroot(fnget_lb2, c(-10, 10), xbxb3gvHbk)
  lb=lbflaginfo$x
  info=1
  
  if (lb < 0 || lb > 1){
    info = 0
  }
  return(c(lb, info))
}
  
