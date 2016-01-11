#' Obtains scaled mean age at death for short growth periods
#'
#' @description Obtains scaled mean age at death assuming a short growth period relative to the life span
#' @family scaled get functions
#' @param pars 4 or 7-vector with parameters: [g lT ha sG] or [g k lT vHb vHp ha SG]
#' @param f optional scalar with scaled reserve density at birth (default f = 1)
#' @param lb optional scalar with scaled length at birth (default: lb is obtained from get_lb)
#' @param lp optional scalar with scaled length at puberty
#' @details Divide the result by the somatic maintenance rate coefficient to arrive at the mean age at death.
#' The variant get_tm_foetus does the same in case of foetal development.
#' If the input parameter vector has only 4 elements (for [g, lT, ha/ kM2, sG]),
#' it skips the calulation of the survival probability at birth and puberty.
#' @return list with  scalar with scaled mean life span (tm),
#' scalar with survival probability at birth (if length p = 7) (Sb),
#' scalar with survival prabability at puberty (if length p = 7) (Sp)
#' and indicator equals 1 if successful convergence, 0 otherwise (info)
#' @export
get_tm_s = function(pars, f = 1, lb = NA, lp = NA){

  with(as.list(pars), {

    if(abs(sG) < 1e-10)
      sG <- 1e-10;

    li <- f - lT;
    hW3 <- ha * f * g/ 6/ li; hW <- hW3^(1/3); # scaled Weibull aging rate
    hG <- sG * f * g * li^2;  hG3 <- hG^3;     # scaled Gompertz aging rate
    tG <- hG/ hW; tG3 <- hG3/ hW3;             # scaled Gompertz aging rate

    info.lp <- 1;
    if(length(pars) >= 7) {
      if(!is.na(lp) && !is.na(lb))
        list[lp, lb, info.lp] <- get_lp(p, f)
      else if(!is.na(lp))
        list[lp, lb, info.lp] <- get_lp(p, f, lb)

      # get scaled age at birth, puberty: tb, tp
      list[tb, lb, info.tb] <- get_tb(p, f, lb);
      irB <- 3 * (1 + f/ g);
      tp <- tb + irB * log((li - lb)/ (li - lp));
      hGtb <- hG * tb;
      Sb <- exp((1 - exp(hGtb) + hGtb + hGtb^2/ 2) * 6/ tG3);
      hGtp <- hG * tp;
      Sp <- exp((1 - exp(hGtp) + hGtp + hGtp^2/ 2) * 6/ tG3);

      if(info.lp == 1 && info.tb == 1)
        info <- 1
      else
        info <- 0
    } else {   # length(p) == 4
      Sb <- NA; Sp <- NA; info <- 1;
    }

    if(abs(sG) < 1e-10) {
      tm <- gamma(4/3)/ hW;
      tm_tail <- 0;
    } else if(hG > 0) {
        tm = 10/ hG; # upper boundary for main integration of S(t)
        tm_tail = expint(exp(tm * hG) * 6/ tG3)/ hG;
        tm = integrate(function(x) fnget_tm_s(x, tG), 0, tm * hW)$value/ hW + tm_tail;
      } else {  # hG < 0
        tm = -1e4/ hG; # upper boundary for main integration of S(t)
        hw = hW * sqrt( - 3/ tG); # scaled hW
        tm_tail = sqrt(pi)* erfc(tm * hw)/ 2/ hw;
        tm = integrate(function(x) fnget_tm_s(x, tG), 0, tm * hW)$value/ hW + tm_tail;
      }

    return(data.frame(tm, Sb, Sp, info))
  })
}

fnget_tm_s = function(t, tG) {
  # called by get_tm_s for life span at short growth periods
  # integrate ageing surv prob over scaled age
  # t: age * hW
  # S: ageing survival prob

  hGt = tG * t; # age * hG
  if(tG > 0)
    S = exp(-(1 + sum(cumprod(hGt * (1/(4:500))))) * t^3)
  else # tG < 0
    S = exp((1 + hGt + hGt^2/ 2  - exp(hGt)) * 6/ tG^3)
}

expint = function(x) {
  integrate(function(t) exp(-t)/t, x, Inf)$value
}
