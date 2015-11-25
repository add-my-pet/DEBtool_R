#' Computes initial scaled reserve
#'
#' @description Obtains scaled length at birth, given the scaled reserve density at birth. Like get_lbarb, but uses a shooting method in 1 variable.
#' @family scaled get functions
#' @param pars 3-vector with parameters: g, k, vbar_H^b
#' @param eb optional scalar with scaled reserve density at birth (default eb = 1)
#' @return scalar with scaled length at birth (lbarb) and indicator equals 1 if successful, 0 otherwise (info)
#' @examples
#' get_lbarb2(c(g = 10, k = 1, vbarHb = 0.01), 1)
#' @export
get_lbarb2 <- function(pars, eb = NA){
  with(as.list(pars), {

    if (is.na(eb))  eb <- 1

    xb <- g/ (eb + g)
    xb3 <- xb^(1/3)

    rootFunction <- function(var) fnget_lbarb2(var, data.frame(xb, xb3, vbarHb, k))
    rFTop <- rootFunction(1/g) # 1/g is the maximum possible for lbarb

    if (rFTop > 0) {
      sol <- uniroot(rootFunction, c(10^(-7), 1/g), f.upper = rFTop, tol = 1e-8)
      lbarb <- sol$root
      info <- 1
    } else {
      lbarb <- NA
      info <- 0
    }

    return(c(lbarb, info))
  })
}
