#' Filters for allowed parameters of standard DEB model without acceleration
#'
#' @description Checks if parameter values are in the allowable part of the parameter space of
#' standard DEB model without acceleration. Meant to be run in the estimation procedure
#' @family filter functions
#' @param par data frame with parameter values
#' @details The flag is an indicator of reason for not passing the filter and it means
#' 0: parameters pass the filter
#' 1: some parameter is negative
#' 2: some kappa is larger than 1
#' 3: growth efficiency is larger than 1
#' 4: maturity levels do not increase during life cycle
#' 5: puberty cannot be reached
#' @return list with filter and flag
#' @examples filter_std(par)
#' @export
filter_std <- function(par){

  with(as.list(par), {

    filter <- 0; flag <- 0; # default setting of filter and flag

    parvec = c(z, kap.X, kap.P, v, kap, p.M, E.G, k.J, E.Hb, E.Hp, kap.R, h.a, s.G, T.A)
    if(any(parvec <= 0) || p.T < 0) {
      flag <- 1
      return(list(filter, flag))
    }

    if(E.Hb >= E.Hp) {  # maturity at birth, puberty
      flag <- 4
      return(list(filter, flag))
    }

    if(f > 1) {
      flag <- 2
      return(list(filter, flag))
    }

    parvec = c(kap, kap.R, kap.X, kap.P)
    if(any(parvec >= 1)) {
      flag <- 2
      return(list(filter, flag))
    }

    # compute and unpack cpar (compound parameters)
    c = parscomp(par);

    if(c$kap.G >= 1) {
      flag <- 3
      return(list(filter, flag))
    }

    if(c$k * c$v.Hp >= f^3) {
      flag <- 5
      return(list(filter, flag))
    }

    if(!reach_birth(c$g, c$k, c$v.Hb, f)) {
      flag <- 6
      return(list(filter, flag))
    }


    filter <- 1

    return(list(filter, flag))

  })
}
