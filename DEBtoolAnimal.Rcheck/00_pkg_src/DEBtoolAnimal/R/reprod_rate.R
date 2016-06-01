#' Gets reproduction rate
#'
#' @description Calculates the reproduction rate in number of eggs per time for an individual of length L
#' and scaled reserve density f.
#' @family scaled get functions
#' @param L n-vector with length
#' @param f scalar with functional response
#' @param pars 9-vector with parameters: kap, kapR, g, kJ, kM, LT, v, UHb, UHp
#' @param Lf optional scalar with length at birth (initial value only)
#' or optional 2-vector with length, L, and scaled functional response f0
#' for a juvenile that is now exposed to f, but previously at another f
#' @return list with n-vector with reproduction rates (R), scalar with scaled initial reserve (UE0),
#' scalar with (volumetric) length at birth (Lb), scalar with (volumetric) length at puberty (Lp) and
#' indicator with 1 for success, 0 otherwise (info)
#' @export
reprod_rate = function(L, f = 1, pars, Lf = NA){
  #  Explanation of variables:
  #  R = kapR * pR/ E0
  #  pR = (1 - kap) pC - kJ * EHp
  #  [pC] = [Em] (v/ L + kM (1 + LT/L)) f g/ (f + g); pC = [pC] L^3
  #  [Em] = {pAm}/ v
  #
  #  remove energies; now in lengths, time only
  #
  #  U0 = E0/{pAm}; UHp = EHp/{pAm}; SC = pC/{pAm}; SR = pR/{pAm}
  #  R = kapR SR/ U0
  #  SR = (1 - kap) SC - kJ * UHp
  #  SC = f (g/ L + (1 + LT/L)/ Lm)/ (f + g); Lm = v/ (kM g)

  with(as.list(pars), {
    Lm <- v/ (kM * g); # maximum length
    k <- kJ/ kM;       # -, maintenance ratio
    VHb <- UHb/ (1 - kap); VHp <- UHp/ (1 - kap);
    vHb <- VHb * g^2 * kM^3/ v^2; vHp <- VHp * g^2 * kM^3/ v^2;

    pars.lb <- c(g = g, k = k, vHb = vHb); # pars for get_lb
    pars.lp <- c(g = g, k = k, lT = LT/ Lm, vHb = vHb, vHp = vHp); # pars for get_tp
    pars.UE0 <- c(VHb = VHb, g = g, kJ = kJ, kM = kM, v = v);      # pars for initial_scaled_reserve
    pars.mat <- c(kap = kap, kapR = kapR, g = g, kJ = kJ, kM = kM, LT = LT, v = v, UHb = UHb, UHp = UHp); # pars for maturity

    if(length(Lf) <= 1) {
      lb0 <- Lf/ Lm; # scaled length at birth
      list[lp, lb, info.lp] <- get_lp(pars.lp, f, lb0);
      Lb <- lb * Lm; Lp <- lp * Lm; # volumetric length at birth, puberty
      if(info.lp != 1) { # return at failure for tp
        cat("lp could not be obtained in reprod_rate \n")
        return(list(R <- L * 0, UE0 <- NA, Lb = Lb, Lp = Lp, info))
      }
    } else { # if length Lb0 = 2
      L0 <- Lf(1); # cm, structural length at time 0
      f0 <- Lf(2); # -, scaled func response before time 0
      list[UH0, info.mat] = maturity(L0, f0, pars.mat);  # d.cm^2, maturity at zero
      if(info.mat != 1) { # return at failure for tp
        cat("maturity could not be obtained in reprod_rate \n")
        return(list(R <- L * 0, UE0 <- NA, Lb = Lb, Lp = Lp, info))
      }
      list[lp, lb, info_lp] <- get_lp(pars.lp, f)
      Lb <- lb * Lm
      UHspan <- seq(UH0, UHp, length=100)
      tLini <- c(0, L0)
      sol <- ode(tLini, UHspan, dget_tL, parms = data.frame(f, g, v, kap, kJ, Lm, LT), method="ode45")
      lp <- min(li - 1e-4, tail(sol[,2], n = 1));

      #list[UH, tL] <- ode45(@dget_tL, [UH0; UHp], [0; L0], [], f, g, v, kap, kJ, Lm, LT);
      Lp <- tL(end,2);  # cm, struc length at puberty after time 0
    }

    list[UE0, Lb, info] <- initial_scaled_reserve(f, pars.UE0, Lb)
    SC <- f * L^3 * (g/ L + (1 + LT/ L)/ Lm)/ (f + g)
    SR <- (1 - kap) * SC - kJ * UHp
    R <- (L >= Lp) * kapR * SR/ UE0  # set reprod rate of juveniles to zero

    return(list(R <- R, UE0 <- UE0, Lb = Lb, Lp = Lp, info))
  })
}



dget_tL <-function (UH, tL, f, g, v, kap, kJ, Lm, LT) {
  # called by cum_reprod
  L <- tL[2]
  r <- v * (f/ L - (1 + LT/ L)/ Lm)/ (f + g)  # 1/d, spec growth rate
  dL <- L * r/ 3  # cm/d, d/dt L
  dUH <- (1 - kap) * L^3 * f * (1/ L - r/ v) - kJ * UH  # cm^2, d/dt UH
  dtL <- c(1, dL)/ dUH  # 1/cm, d/dUH L
  return(dtL)
}

