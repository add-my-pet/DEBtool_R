#' Warns of unreasonable parameters for the standard DEB model without acceleration
#'
#' @description Checks if parameter values are in the reasonable part of the parameter space
#' of standard DEB model without acceleration, produces warnings. Meant to be run after the estimation procedure
#' @family filter functions
#' @param par data frame with parameter values
#' @examples warning_std(par)
#' @export
warning_std <- function(par){

  compPar <- parscomp(par);

  with(as.list(par), with(as.list(compPar), {

    if(exists("kap.P"))
      if(kap.X + kap.P >= 1)
        cat("kap_X + kap_P > 1, which violates energy conservation. \n");

    if(kap.G >= mu.V / mu.E)
      cat("kap_G >= mu_V / mu_E, which is not allowed if CO2 production occurs in association with growth. \n");

    if(kap.X >= mu.E / mu.X)
      cat("kap_X > mu_X / mu_E, which is not allowed if CO2 production occurs in association with assimilation. \n");

    if(exists("kap.P"))
      if(kap.X * mu.X / mu.E + kap.P * mu.X / mu.P >= 1)
        cat("kap_X * mu_X / mu_E + kap_P * mu_X / mu_P > 1, which is not allowed if CO2 production occurs in association with assimilation. \n");

    pars.tp <- c(g, k, l.T, v.Hb, v.Hp)
    t.p <- get_tp(pars.tp, 1)
    pars.tm <- c(g, l.T, h.a/ k.M^2, s.G)
    t.m <- get_tm_s(pars.tm, 1)

    if(t.m < t.p)
      cat("Ageing is too fast for the organism to be able to reproduce if death occurs at (mean) life span. \n")

  }))
}

