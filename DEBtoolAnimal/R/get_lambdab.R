#' Computes initial scaled reserve
#'
#' @description Obtains scaled length at birth, given the scaled reserve density at birth.
#' @family scaled get functions
#' @param p 3-vector with parameters: g, k, vv_H^b (see below)
#' @param eb optional scalar with scaled reserve density at birth (default eb = 1)
#' @param lambdab0 optional scalar with initial estimate for scaled length at birth (default lambdab0: lambdab for k = 1)
#' @return scalar with scaled length at birth (lambdab) and indicator equals 1 if successful, 0 otherwise (info)
#' @examples
#' beta0(0.1, 0.2)
#' @export
get_lambdab = function(p, eb, lambdab0=NA){

  ## Remarks
  # The theory behind get_lb, get_tb and get_ue0 is discussed in
  #    <http://www.bio.vu.nl/thb/research/bib/Kooy2009b.html Kooy2009b>.
  # Solves y(x_b) = y_b  for lb with explicit solution for y(x)
  #   y(x) = x e_H/(1-kap) = x g u_H/ l^3
  #   and y_b = x_b g u_H^b/ ((1-kap)l_b^3)
  #   d/dx y = r(x) - y s(x);
  #   with solution y(x) = v(x) \int r(x)/ v(x) dx
  #   and v(x) = exp(- \int s(x) dx).
  # A Newton Raphson scheme is used with Euler integration, starting from an optional initial value.
  # Replacement of Euler integration by ode23: <get_lb1.html *get_lb1*>,
  #  but that function is much lower.
  # Shooting method: <get_lb2.html *get_lb2*>.
  # Bisection method (via fzero): <get_lb3.html *get_lb3*>.
  # In case of no convergence, <get_lb2.html *get_lb2*> is run automatically as backup.
  # Consider the application of <get_lb_foetus.html *get_lb_foetus*> for an alternative initial value.

  #  unpack p
  g <- p[1]   # g = [E_G] * v/ kap * {p_Am}, energy investment ratio
  k <- p[2]   # k = k_J/ k_M, ratio of maturity and somatic maintenance rate coeff
  vvHb <- p[3] # vv_H^b = v_H^b/ g^3; v_H^b = U_H^b g^2 kM^3/ (1 - kap) v^2; U_H^b = M_H^b/ {J_EAm}

  lambdab0 <- Re(lambdab0)

  info = 1
  if (!exists('lb0')){
    lb0 = NA
  }

  if (k == 1){
    lambdab = as.complex(vvHb)^(1/ 3) # exact solution for k = 1
    info = 1
  }

  if (is.na(lb0)){
    lambdab = as.complex(vvHb)^(1/3) # exact solution for k = 1
  } else {lambdab = lambdab0}

  if (!exists('eb')){
    eb = 1
  } else if (is.na(eb)){
    eb = 1
  }

  n = 1000 + round(1000 * max(0, k - 1))
  xb = g/ (g + eb)
  xb3 = xb ^ (1/3)
  x = seq(1e-6, xb, length=n)
  dx = xb/ n
  x3 = x ^ (1/3)

  b = beta0(x, xb)/ 3
  t0 = xb * vvHb
  i = 0
  norm = 1 # make sure that we start Newton Raphson procedure
  ni = 100 # max number of iterations

  print(lambdab)

  while (i < ni  && norm > 1e-8 && !is.na(norm)){
    print(x3, xb3)
    lambda = x3 / (xb3/ lambdab - b)
    s = (k - x) / (1 - x) * lambda / x
    v = exp( - dx * cumsum(s))
    vb = v[n]
    rho = 1 + lambda
    rhov = rho / v
    t = t0/ lambdab^3/ vb - dx * sum(rhov)
    dl = xb3/ lambdab^2 * lambda ^ 2 / x3
    dlnl = dl / lambdab
    dv = v * exp( - dx * cumsum(s * dlnl))
    dvb = dv[n]
    dlnv = dv / v
    dlnvb = dlnv[n]
    drho = dl
    dlnrho = drho / rho
    dt = - t0/ lambdab^3/ vb * (3/ lambdab + dlnvb) - dx * sum((dlnrho - dlnv) * rhov)
    # [i lb t dt] # print progress
    lambdab = Re(lambdab - t/ dt) # Newton Raphson step
    norm = Re( t^2 )
    i = i + 1
  }

#  if (i == ni || lb < 0 || lb > 1 || is.na(norm)) { # no convergence
#    # try to recover with a shooting method
#    if (is.na(lb0)){
#      lbinfo = get_lb2(p, eb)
#      lb=lbinfo[1]
#      info=lbinfo[2]
#    } else if (lb0 < 1 && lb0 > 0){
#      lbinfo = get_lb2(p, eb, lb0)
#      lb=lbinfo[1]
#      info=lbinfo[2]
#    } else {
#      lbinfo = get_lb2(p, eb)
#      lb=lbinfo[1]
#      info=lbinfo[2]
#    }
#  }

  if (info == 0){
    print('warning get_lb: no convergence of l_b')
  }

  return(c(lambdab, info))

}
