#' Gets initial scaled reserve
#'
#' @description Obtains scaled length at birth, given the scaled reserve density at birth.
#' @family scaled get functions
#' @param f n-vector with scaled functional responses
#' @param pars 3-vector with parameters: g, k, vv_H^b (see below)
#' @param lb0 optional n-vector with lengths at birth
#' @return n-vector with initial scaled reserve: M_E^0/ {J_EAm} (U0), n-vector with length at birth 8Lb) and n-vector with 1's if successful, 0's otherwise (info)
#' @examples
#' get_lambdab(c(g = 10, k = 1, vvHb = 0.0005), 1)
initial_scaled_reserve <- function(f, p, Lb0){

  ## Description
  # Gets initial scaled reserve
  #
  # Input
  #
  # * f: n-vector with scaled functional responses
  # * p: 5-vector with parameters: VHb, g, kJ, kM, v
  # * Lb0: optional n-vector with lengths at birth
  #
  # Output
  #
  # * U0: n-vector with initial scaled reserve: M_E^0/ {J_EAm}
  # * Lb: n-vector with length at birth
  # * info: n-vector with 1's if successful, 0's otherwise

  ## Remarks
  # Like <get_ue0.html *get_ue0*>, but allows for vector arguments and
  # input and output is not downscaled to dimensionless quantities,

  ## Example of use
  # p = [.8 .42 1.7 1.7 3.24 .012]; initial_scaled_reserve(1,p)

  #  unpack parameters
  VHb = p[1] # d mm^2, scaled maturity at birth: M_H^b/[[1-kap]{J_EAm}]
  g   = p[2] # -, energy investment ratio
  kJ  = p[3] # 1/d, maturity maintenance rate coefficient
  kM  = p[4] # 1/d, somatic maintenance rate coefficient
  v   = p[5] # mm/d, energy conductance

  # if kJ = kM: VHb = g * Lb^3/ v;

  nf = length(f)
  U0 = rep(0, length=nf)
  Lb = rep(0, length=nf)
  info = rep(0, length=nf)
  q = c(g, kJ/ kM, VHb * g^2 * kM^3/ v^2)
  if (exists('Lb0')){
    lb0 = rep(1, length=nf) * Lb0 * kM * g/ v
  } else {
    lb0 = ones(nf,1) * get_lb(q,f[1])[1] # initial estimate for scaled length
  }
  for (i in c(1:nf)){
    lbinfo = get_lb(q, f[i], lb0[i])
    lb=lbinfo[1]
    info[i]=lbinfo[2]
    ## try get_lb1 or get_lb2 for higher accuracy
    Lb[i] = lb * v/ kM/ g;
    uE0lbinfo = get_ue0(q, f[i], lb)
    uE0=uE0lbinfo[1]
    U0[i] = uE0 * v^2/ g^2/ kM^3
  }

  return(c(U0,Lb, info))
}
