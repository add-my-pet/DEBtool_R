#' Gets initial scaled reserve
#'
#' @description Gets initial scaled reserve.
#' @family scaled get functions
#' @param f n-vector with scaled functional responses
#' @param pars 5-vector with parameters: VHb, g, kJ, kM, v
#' @param Lb0 optional n-vector with lengths at birth
#' @return n-vector with initial scaled reserve: M_E^0/ {J_EAm} (U0), n-vector with length at birth (Lb) and n-vector with 1's if successful, 0's otherwise (info)
#' @examples
#' initial_scaled_reserve(f = c(1, 0.9), pars = c(VHb = .8, g = .42, kJ = 1.7, kM = 1.7, v = 3.24))
#' @export
initial_scaled_reserve <- function(f, pars, Lb0 = NA){

  with(as.list(pars), {

  # if kJ = kM: VHb = g * Lb^3/ v;

    nf <- length(f)
    U0 <- rep(0, length=nf)
    Lb <- rep(0, length=nf)
    info <- rep(0, length=nf)
    q = c(g = g, k = kJ/ kM, vHb = VHb * g^2 * kM^3/ v^2)
    if (length(Lb0) == 1 && is.na(Lb0)) {
      lb0 <- rep(1, length=nf) * get_lb(p = q, eb = f[1]) # initial estimate for scaled length
    } else {
      lb0 <- rep(1, length=nf) * Lb0 * kM * g/ v
    }
    for (i in c(1:nf)){
      lbinfo <- get_lb(q, f[i], lb0[i])
      lb <- lbinfo[1]
      info[i] <- lbinfo[2]
      ## try get_lb1 or get_lb2 for higher accuracy
      Lb[i] <- lb * v/ kM/ g;
      uE0lbinfo <- get_ue0(q, f[i], lb)
      uE0 <- uE0lbinfo[[1]]
      U0[i] <- uE0 * v^2/ g^2/ kM^3
    }

    return(data.frame(U0, Lb, info))
  })
}
