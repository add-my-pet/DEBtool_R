#' Computes initial scaled reserve
#'
#' @description Obtains the initial scaled reserve given the scaled reserve density at birth.
#' Function get_ue0 does so for eggs, get_ue0_foetus for foetuses.
#' Specification of length at birth as third input by-passes its computation,
#' so if you want to specify an initial value for this quantity, you should use get_lb directly.
#' @family get functions
#' @param pars 1 or 3 -vector with parameters g, k_J/ k_M, v_H^b, see get_lb
#' @param eb optional scalar with scaled reserbe density at birth (default: eb = 1)
#' @param lb0 optional scalar with scaled length at birth (default: lb is optained from get_lb)
#' @return uE0 scalar with scaled reserve at t=0: $U_E^0 g^2 k_M^3/ v^2$ with $U_E^0 = M_E^0/ \{J_{EAm}\}$, lb scalar with scaled length at birth and info indicator equals 1 if successful, 0 otherwise
#' @examples
#' get_ue0(pars = c(0.42, 1, 0.066), eb = 1, lb0 = 0.4042)
#' @export
get_ue0 <- function(pars, eb = 1, lb0 = NA){
  with(as.list(pars), {
    if (is.na(lb0)){
      if (length(p) < 3) {
        print("not enough input parameters, see get_lb")
        return(c(uE0 = NA, lb = NA, info = 0))
      }
      lbinfo = get_lb(p = pars, eb = eb )
      lb = lbinfo[1]
      info = lbinfo[2]
    } else {
      lb = lb0
      info = 1
    }

    xb = g/ (eb + g)
    uE0 = (3 * g/ (3 * g * xb^(1/ 3)/ lb - beta0(0, xb)))^3

    return(c(uE0, lb, info))
  })
}
