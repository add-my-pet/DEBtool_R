#' Computes f using the ode solver for delta(x), for finding lambdab
#'
#' @description Computes f using the ode solver for delta(x), for finding lambdab.
#' @family scaled get functions
#' @param lambdab scalar with scaled length at birth (lambdab = lb/ g)
#' @param pars data.frame with lambdab, xb, xb3 (xb^1/3), k
#' @return scalar with function f which when zero indicates lambdab
#' @examples
#' fnget_lambdab2(0.03, c(xb = 10/11, xb3 = (10/11)^(1/3), vvHb = 0.001, k = 1))
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
