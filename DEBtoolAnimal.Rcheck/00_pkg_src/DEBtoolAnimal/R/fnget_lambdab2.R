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
fnget_lambdab2 <- function(lambdab, pars){
  # f = y(x_b) - y_b = 0; x = g/ (e + g); x_b = g/ (e_b + g)
  # y(x) = x e_H = x g u_H/ l^3 and y_b = x_b g u_H^b/ l_b^3
  with(as.list(pars), {
    tspan <- seq(1e-10, xb, length=100)
    deltaini <- c(delta = 0)
    xdelta <- ode(deltaini, tspan, dget_lambdab2, parms = data.frame(lambdab, xb, xb3, k), method="ode23")
    xdelta[nrow(xdelta),2] - xb * vvHb/ lambdab^3  # f
  })
}
