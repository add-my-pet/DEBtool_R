#' Computes initial scaled reserve density at birth
#'
#' @description Obtains the initial scaled reserve given the scaled reserve density at birth.
#' Function get_ue0 does so for eggs, get_ue0_foetus for foetuses.
#' Specification of length at birth as third input by-passes its computation,
#' so if you want to specify an initial value for this quantity, you should use get_lb directly.
#' @family get functions
#' @param eb: optional scalar with scaled reserbe density at birth
#   (default: eb = 1)
#' @param x1 scalar with upper boundary for integration
#' @return scalar with particular incomple beta function
#' @examples
#' get_nue0(g = 10, lambdab = 0.01)
#' get_nue0(g = 10, k = 0.7, vvHb = 5e-4)
#' @export

## get_ue0
# gets initial scaled reserve

##
#  function [uE0, lb, info] =
get_nue0 <- function(g = NA, k = NA, vvHb = NA, eb = 1, lambdab = NA){
  # Output
  #
  # * uE0: scaled with scaled reserve at t=0: $U_E^0 g^2 k_M^3/ v^2$
    #   with $U_E^0 = M_E^0/ \{J_{EAm}\} or E^0/ \{p_{Am}\}$
    # * lb: scalar with scaled length at birth
  # * info: indicator equals 1 if successful, 0 otherwise

  if (is.na(g))
    stop("input value for g needed")

  xb = g/ (eb + g)

  if (!is.na(lambdab)) {
    info = 1
  } else if(!is.na(k) && !is.na(vvHb)) {
    #list[lambdab, info] <- get_lambdab(c(g = g, k = k, vvHb = vvHb), eb)
    sol <- get_lambdab(c(g = g, k = k, vvHb = vvHb), eb)
    lambdab <- sol[1]
    info <- sol[2]

  } else
    stop("input value for lambdab or (k and vvHb) needed")

  nuE0 = (3/ (3 * xb^(1/ 3)/ lambdab - beta0(0, xb)))^3

  return(c(nuE0, lambdab, info))
}
