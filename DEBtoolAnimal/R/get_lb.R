#' Computes scaled length at birth
#'
#' @description Obtains scaled length at birth, given the scaled reserve density at birth.
#' @family scaled get functions
#' @param pars 3-vector with parameters: g, k, v_H^b
#' @param eb optional scalar with scaled reserve density at birth (default eb = 1)
#' @param lbarb0 optional scalar with initial estimate for scaled length at birth (default lb0: lb for k = 1)
#' @return scalar with scaled length at birth (lb) and indicator equals 1 if successful, 0 otherwise (info)
#' @examples
#' get_lb(c(g = 10, k = 1, vHb = 0.5), 1)
#' @export
get_lb = function(pars, eb = 1, lb0 = as.numeric(pars[3]^(1/3))){
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

  with(as.list(pars), {
    info <- 1

    if (k == 1){
      lb <- vHb^(1/3) # exact solution for k = 1
      return(c(lb, info))
    }

    lb <- lb0

    n <- 1000 + round(1000 * max(0, k - 1))
    xb <- g/ (g + eb)
    xb3 <- xb ^ (1/3)
    x <- seq(1e-6, xb, length=n)
    dx <- xb/ n
    x3 <- x ^ (1/3)

    b = beta0(x, xb)/ (3 * g)
    t0 = xb * g * vHb
    i = 0
    norm = 1 # make sure that we start Newton Raphson procedure
    ni = 100 # max number of iterations

    while (i < ni  && norm > 1e-8 && !is.na(norm)){
      l = x3 / (xb3/ lb - b)
      s = (k - x) / (1 - x) * l/ g / x
      v = exp( - dx * cumsum(s))
      vb = v[n]
      r = (g + l)
      rv = r / v
      t = t0/ lb^3/ vb - dx * sum(rv)
      dl = xb3/ lb^2 * l ^ 2 / x3
      dlnl = dl / l
      dv = v * exp( - dx * cumsum(s * dlnl))
      dvb = dv[n]
      dlnv = dv / v
      dlnvb = dlnv[n]
      dr = dl
      dlnr = dr / r
      dt = - t0/ lb^3/ vb * (3/ lb + dlnvb) - dx * sum((dlnr - dlnv) * rv)
      lb = Re(lb - t/ dt) # Newton Raphson step
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
      print("warning get_lb: no convergence of l_b")
    }

    return(data.frame(lb, info))
  })
}
