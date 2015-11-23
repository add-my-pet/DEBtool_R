#' Computes derivative d delta/dx
#'
#' @description Obtains the derivative d delta/dx from lambdab, xb and k.
#' @family scaled get functions
#' @param x scalar x = g/(g + e)
#' @param delta scalar delta = x e_H/ (1 - kap)g
#' @param pars data.frame with lambdab, xb, xb3 (xb^1/3), k
#' @return scalar with derivative value d delta/ dx
#' @examples
#' dget_lambdab2(10^(-6), 0, c(lambdab = 0.003, xb = 10/11, xb3 = (10/11)^(1/3), k = 1))
#' @export
dget_lambdab2=function(x, delta, pars){
  # delta = y/ g; y = x e_H/ (1 - kap); x = g/(g + e)
  # (x,delta): (0,0) -> (xb, xb eHb/ (1 - kap)/ g)
  with(as.list(pars), {
    x3 <- x^(1/ 3)
    lambda <- x3/ (xb3/ lambdab - beta0(x, xb)/ 3)
    ddelta <- 1 + lambda - delta * (k - x)/ (1 - x) * lambda/ x
    return(list(ddelta))
  })
}
