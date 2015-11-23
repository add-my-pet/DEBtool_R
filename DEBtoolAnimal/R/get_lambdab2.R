#' Computes initial scaled reserve
#'
#' @description Obtains scaled length at birth, given the scaled reserve density at birth. Like get_lambdab, but uses a shooting method in 1 variable.
#' @family scaled get functions
#' @param pars 3-vector with parameters: g, k, vv_H^b (see below)
#' @param eb optional scalar with scaled reserve density at birth (default eb = 1)
#' @return scalar with scaled length at birth (lambdab) and indicator equals 1 if successful, 0 otherwise (info)
#' @examples
#' get_lambdab(c(g = 10, k = 1, vvHb = 0.01), 1)
#' @export
get_lambdab2 <- function(pars, eb = NA){
  with(as.list(pars), {

    if (!exists('eb') || is.na(eb))  eb <- 1

    xb <- g/ (eb + g)
    xb3 <- xb^(1/3)

    rootFunction <- function(var) fnget_lambdab2(var, data.frame(xb, xb3, vvHb, k))
    rFTop <- rootFunction(1/g) # 1/g is the maximum possible for lambdab

    if (rFTop > 0) {
      sol <- uniroot(rootFunction, c(10^(-7), 1/g), f.upper = rFTop, tol = 1e-8)
      lambdab <- sol$root
      info <- 1
    } else {
      lambdab <- NA
      info <- 0
    }

    return(c(lambdab, info))
  })
}
