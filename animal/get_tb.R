## get_tb
#

##
get_tb= function(p, eb, lb=NA){
#  created at 2007/07/27 by Bas Kooijman; modified 2014/03/17, 2015/01/18

## Syntax
# [tb, lb, info] = <../get_tb.m *get_tb*> (p, eb, lb)

## Description
# Obtains scaled age at birth, given the scaled reserve density at birth. 
# Divide the result by the somatic maintenance rate coefficient to arrive at age at birth. 
#
# Input
#
# * p: 1 or 3-vector with parameters g, k_J/ k_M, v_H^b
#
#     Last 2 values are optional in invoke call to get_lb
#
# * eb: optional scalar with scaled reserve density at birth (default eb = 1)
# * lb: optional scalar with scaled length at birth (default: lb is obtained from get_lb)
#  
# Output
#
# * tb: scaled age at birth \tau_b = a_b k_M
# * lb: scalar with scaled length at birth: L_b/ L_m
# * info: indicator equals 1 if successful, 0 otherwise

## Remarks
# See also <get_tb1.html *get_tb1*> for backward integration over maturity and
# <get_tb_foetus.html *get_tb_foetus*> for foetal development

## Example of use
# get_tb([.1;.5;.03])
# See also <../ mydata_ue0.m *mydata_ue0*>
  
  if (!exists('eb')) {
    eb = 1                   # maximum value as juvenile
  }
  
  info = 1
  
  if (!exists('lb')) {
    if (length(p) < 3){
      print('not enough input parameters, see get_lb \n')
      tb = NA
    }
    lbinfo = get_lb(p, eb)
    lb=lbinfo[1]
    info=lbinfo[2]
  }
  if (is.na(lb)){
    lbinfo = get_lb(p, eb)
    lb=lbinfo[1]
    info=lbinfo[2]
  }
  
  # unpack p
  g = p[1]  # energy investment ratio
  
  xb = g/ (eb + g)  # f = e_b 
  ab = 3 * g * xb^(1/ 3)/ lb  # \alpha_b
  library(pracma)
  abxb=c(ab,xb)
  tb = 3 * quad(dget_tb, xa= 1e-15, xb=xb, tol = 1.0e-12, trace = FALSE, abxb)
  return(c(tb, lb, info))
}
  
  # subfunction
  
dget_tb= function(x, abxb){
  # called by get_tb
  ab=abxb[1]
  xb=abxb[2]
  f = x ^ (-2/ 3) / (1 - x) / (ab - beta0(x, xb))
  return(f)
}

