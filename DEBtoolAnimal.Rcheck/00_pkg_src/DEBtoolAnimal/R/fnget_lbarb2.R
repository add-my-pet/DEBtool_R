#' Computes f using the ode solver for delta(x), for finding lbarb
#'
#' @description Computes f using the ode solver for delta(x), for finding lbarb.
#' @family scaled get functions
#' @param lbarb scalar with scaled length at birth (lbarb = lb/ g)
#' @param pars data.frame with lbarb, xb, xb3 (xb^1/3), k
#' @return scalar with function f which when zero indicates lbarb
#' @examples
#' fnget_lbarb2(0.03, c(xb = 10/11, xb3 = (10/11)^(1/3), vbarHb = 0.001, k = 1))
#' @export
fnget_lbarb2 <- function(lbarb, pars){
  # f = y(x_b) - y_b = 0; x = g/ (e + g); x_b = g/ (e_b + g)
  # y(x) = x e_H = x g u_H/ l^3 and y_b = x_b g u_H^b/ l_b^3
  with(as.list(pars), {
    tspan <- seq(1e-10, xb, length=100)
    deltaini <- c(delta = 0)
    xdelta <- ode(deltaini, tspan, dget_lbarb2, parms = data.frame(lbarb, xb, xb3, k), method="ode23")
    xdelta[nrow(xdelta),2] - xb * vbarHb/ lbarb^3  # f
  })
}
