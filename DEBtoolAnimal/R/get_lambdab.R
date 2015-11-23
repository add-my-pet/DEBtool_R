#' Computes initial scaled reserve
#'
#' @description Obtains scaled length at birth, given the scaled reserve density at birth.
#' @family scaled get functions
#' @param pars 3-vector with parameters: g, k, vv_H^b (see below)
#' @param eb optional scalar with scaled reserve density at birth (default eb = 1)
#' @param lambdab0 optional scalar with initial estimate for scaled length at birth (default lambdab0: lambdab for k = 1)
#' @return scalar with scaled length at birth (lambdab) and indicator equals 1 if successful, 0 otherwise (info)
#' @examples
#' get_lambdab(c(g = 10, k = 1, vvHb = 0.0005), 1)
#' @export
get_lambdab <- function(pars, eb = 1, lambdab0 = NA){

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
    if (!exists('lambdab0'))
      lambdab0 <- NA

    if (k == 1) {
      lambdab <- vvHb^(1/ 3) # exact solution for k = 1
      return(c(lambdab, info))
    }

    if (is.na(lambdab0))
      lambdab <- vvHb^(1/3) # exact solution for k = 1
    else
      lambdab <- lambdab0

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

    while (i < ni && norm > 1e-8 && !is.na(norm)){
      lambda = x3 / (xb3/ lambdab - b)
      s = (k - x) / (1 - x) * lambda / x
      v = exp( - dx * cumsum(s))
      vb = v[n]
      rho = 1 + lambda
      rhov = rho / v
      t = t0/ lambdab^3/ vb - dx * sum(rhov)
      dl = xb3/ lambdab^2 * lambda ^ 2 / x3;     dlnl = dl / lambda
      dv = v * exp( - dx * cumsum(s * dlnl));    dlnv = dv / v
      dvb = dv[n];                               dlnvb = dlnv[n]
      drho = dl;                                 dlnrho = drho / rho
      dt = - t0/ lambdab^3/ vb * (3/ lambdab + dlnvb) - dx * sum((dlnrho - dlnv) * rhov)
      lambdab = Re(lambdab - t/ dt) # Newton Raphson step
      norm = Re( t^2 )
      i = i + 1
    }

    if (i == ni || lambdab < 0 || lambdab > 1/g || is.na(norm)) { # no convergence
      # try to recover with a shooting method
      lambdabinfo <- get_lambdab2(c(g = g, k = k, vvHb = vvHb), eb)
      info    <- lambdabinfo[2]
      if (info == 1)
        lambdab <- lambdabinfo[1]
      else
        print("warning get_lb: no convergence of l_b")
    }

    return(c(lambdab, info))
  })
}
