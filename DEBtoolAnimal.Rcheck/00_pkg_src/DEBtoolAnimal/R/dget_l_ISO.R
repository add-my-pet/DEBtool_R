#' Computes derivative d l/d vH
#'
#' @description Obtains the derivative d l/d vH from g, k, lT, f, sM.
#' @family scaled dget functions
#' @param vH scaled maturity volume
#' @param l scaled length
#' @param pars data.frame with g, k, lT, f, sM
#' @return scalar with derivative value d l/ d vH
#' @export
dget_l_ISO <- function(vH, l, pars){
  with(as.list(pars), {
    r <- g * (f * sM - lT * sM - l)/ l/ (f + g); # specific growth rate
    dl <- l * r/ 3;                              # d/dt l
    dvH <- f * l^2 * (sM - l * r/ g) - k * vH;   # d/dt vH
    return(list(dl = dl/ dvH))
  })
}
