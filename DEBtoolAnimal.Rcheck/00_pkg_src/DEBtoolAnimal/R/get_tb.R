#' Gets scaled age at birth
#'
#' @description Obtains scaled age at birth, given the scaled reserve density at birth.
#' @family scaled get functions
#' @param pars 3-vector with parameters: g, k, v_H^b
#' @param eb optional scalar with scaled reserve density at birth (default eb = 1)
#' @param lb0 optional scalar with initial estimate for scaled length at birth (default lb0: lb for k = 1)
#' @details Multiply the result with the somatic maintenance rate coefficient to arrive at age at puberty.
#' @return list with scaled age at birth tau_b = a_b k_M (tb), scaled length at birth (lb)
#' and indicator equals 1 if successful convergence, 0 otherwise (info)
#' @examples
#' get_tb(c(g = 10, k = 1, vHb = 0.5), 1)
#' @export
get_tb = function(pars, eb = 1, lb = NA){

  with(as.list(pars), {

    if(is.na(lb)) {
      if(length(pars) < 3) {
        cat("not enough input parameters, see get_lb \n");
        return(list(NA, lb, info))
      }
      list[lb, info] <- get_lb(pars, eb);
    }

    xb <- g/ (eb + g); # f = e_b
    ab <- 3 * g * xb^(1/ 3)/ lb   # \alpha_b
    tb <- 3 * integrate(function(x) dget_tb(x, ab, xb), 1e-15, xb)$value

    return(data.frame(tb, lb, info <- 1))
  })
}

dget_tb = function(x, ab, xb)
  # called by get_tb
  return(x ^ (-2/ 3) / (1 - x) / (ab - beta0(x, xb)))
