#' Computes scaled length at puberty
#'
#' @description Obtains scaled length at pubertyat constant food density.
#' @family scaled get functions
#' @param pars 5-vector with parameters: g, k, l_T, v_H^b, v_H^p
#' @param eb optional scalar with scaled reserve density at birth (default eb = 1)
#' @param lb0 optional scalar with initial estimate for scaled length at birth (default lb0: lb for k = 1)
#' @details If scaled length at birth (second input) is not specified, it is computed (using automatic initial estimate).
#' If it is specified, however, is it just copied to the (second) output. Food density is assumed to be constant.
#' @return scaled length at puberty (lp), scaler length at birth (lb)
#' and indicator equals 1 if successful convergence, 0 otherwise (info)
#' @examples
#' get_lp(c(g = 10, k = 1, lT = 0, vHb = 0.5, vHp = 10), 1)
#' @export
get_lp = function(pars, f = 1, lb0 = NA){

  with(as.list(pars), {
    pars.lb <- c(g = g, k = k, vHb = vHb);
    lp <- list(); info <- 0;

    li <- f - lT; # -, scaled ultimate length

    if(is.na(lb0))
      list[lb, info] <- get_lb(pars.lb, f)
    else if(length(lb0) < 2) {
      info <- 1
      lb <- lb0
    } else {  # for a juvenile of length l and maturity vH
      l <- lb0[1]; vH <- lb0[2]; # juvenile now exposed to
      list[lb, info] <- get_lb(pars.lb, f)
    }

    # d/d tau vH = b2 l^2 + b3 l^3 - k vH
    b3 <- 1 / (1 + g / f);
    b2 <- f - b3 * li;
    # vH(t) = - a0 - a1 exp(- rB t) - a2 exp(- 2 rB t) - a3 exp(- 3 rB t) +
    #         + (vHb + a0 + a1 + a2 + a3) exp(-kt)
    a0 <- - (b2 + b3 * li) * li^2/ k; # see get_lp1
    if(vHp > -a0) {     # vH can only reach -a0
      cat("Warning in get_lp: maturity at puberty cannot be reached \n");
      return(data.frame(lp, lb, info))
    }

    if(k == 1)
      lp <- vHp^(1/3)
    else if(length(lb0) < 2)
      if (f * lb^2 * (g + lb) < vHb * k * (g + f))
        cat("Warning in get_lp: maturity does not increase at birth \n")
      else {
        vHspan <- seq(vHb, vHp, length=100)
        lini <- c(l = lb)
        sol <- ode(lini, vHspan, dget_l_ISO, parms = data.frame(k, lT, g, f, sM = 1), method="ode45")
        lp <- min(li - 1e-4, tail(sol[,2], n = 1));
      }
    else # for a juvenile of length l and maturity vH
      if(f * l^2 * (g + l) < vH * k * (g + f))
        cat("Warning in get_lp: maturity does not increase initially \n")
      else if(vH + 1e-4 < vHp) {
        vHspan <- seq(vH, vHp, length=100)
        lini <- c(l = l)
        sol <- ode(lini, vHspan, dget_l_ISO, parms = data.frame(k, lT, g, f, sM = 1), method="ode45")
        lp <- tail(sol[,2], n = 1)
      } else {
        lp <- l;
        cat("Warning in get_lp: initial maturity exceeds puberty threshold \n")
      }

    return(data.frame(lp, lb, info))

  })
}
