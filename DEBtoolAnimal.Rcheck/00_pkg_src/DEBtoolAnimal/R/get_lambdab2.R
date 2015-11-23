#' Computes initial scaled reserve
#'
#' @description Obtains scaled length at birth, given the scaled reserve density at birth.
#' @family scaled get functions
#' @param p 3-vector with parameters: g, k, vv_H^b (see below)
#' @param eb optional scalar with scaled reserve density at birth (default eb = 1)
#' @param lambdab0 optional scalar with initial estimate for scaled length at birth (default lambdab0: lambdab for k = 1)
#' @return scalar with scaled length at birth (lambdab) and indicator equals 1 if successful, 0 otherwise (info)
#' @examples
#' get_lambdab(c(10, 1, 0.01), 1, 0.1)
#' @export
get_lambdab2= function(p, eb, lambdab0=NA){
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
  vvHb = p[3] # v_H^b = U_H^b g^2 kM^3/ (1 - kap) v^2; U_H^b = M_H^b/ {J_EAm}

  info = 1

  if (is.na(lambdab0)) {
    lambdab = as.complex(vvHb)^(1/ 3) # exact solution for k = 1
  } else {
    lambdab = lambdab0
  }
  if (!exists('eb')){
    eb = 1
  } else if (is.na(eb)){
    eb = 1
  }


  xb = g/ (eb + g)
  xb3 = xb^(1/3)

  xbxb3vvHbk=c(xb, xb3, vvHb, k)
  lambdab=Re(lambdab)

  rootFunction <- function(var) fnget_lambdab2(var, data.frame(xb, xb3, vvHb, k))
  rFTop <- rootFunction(1/g)

  if (rFTop > 0) {
    #lambdabflaginfo = fzero(fnget_lambdab2, lambdab, maxiter = 100, tol = 10e-8, data.frame(xb, xb3, vvHb, k))
    sol = uniroot(rootFunction, c(10^(-7), 1/g), f.upper = rFTop, tol = 1e-4)
    lambdab = sol$root
    info = 1
  } else {
    info = 0
  }

  return(c(lambdab, info))
}

