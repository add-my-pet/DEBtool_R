#' Obtains predictions, using parameters and data
#'
#' @description Obtains predictions, using parameters and data
#' @family add-my-pet template functions
#' @param par data frame with parameter values
#' @param data data frame with data values
#' @param auxData data frame with auxiliary data values
#' @return list with prdData (data frame with values of predictions) and info (indicator for customized filters)
#' @examples predict_my_pet(par, data, auxData)
#' @export
predict_my_pet <- function(par, data, auxData = list()){

  compPar <- parscomp(par);
  with(as.list(par), with(as.list(compPar), with(as.list(data), with(as.list(auxData), {

    # customized filters for allowable parameters of the standard DEB model (std)
    # for other models consult the appropriate filter function.
    filterChecks <- k * v.Hp >= f.tL^3 ||   # constraint required for reaching puberty with f.tL
      !reach_birth(g, k, v.Hb, f.tL);       # constraint required for reaching birth with f.tL

    if(filterChecks)
      return(prdData <- list(), info <- 0)

    # compute temperature correction factors
    TC.ab <- tempcorr(temp$ab, T.ref, T.A);
    TC.ap <- tempcorr(temp$ap, T.ref, T.A);
    TC.am <- tempcorr(temp$am, T.ref, T.A);
    TC.Ri <- tempcorr(temp$Ri, T.ref, T.A);
    TC.tL <- tempcorr(temp$tL, T.ref, T.A);

    # uncomment if you need this for computing moles of a gas to a volume of gas
    # - else feel free to delete  these lines
    # molar volume of gas at 1 bar and 20 C is 24.4 L/mol
    # T = C2K(20); % K, temp of measurement equipment- apperently this is
    # always the standard unless explicitely stated otherwise in a paper (pers. comm. Mike Kearney).
    # X.gas = T.ref/ T/ 24.4;  # M,  mol of gas per litre at T.ref and 1 bar;

    # zero-variate data

    # life cycle
    pars.tp <- c(g = g, k = k, lT = l.T, vHb = v.Hb, vHp = v.Hp);               # compose parameter vector
    list[t.p, t.b, l.p, l.b, info] <- get_tp(pars.tp, f); # -, scaled times & lengths at f

    # birth
    L.b <- L.m * l.b;                  # cm, structural length at birth at f
    Lw.b <- L.b/ del.M;                # cm, physical length at birth at f
    Wd.b <- L.b^3 * d.V * (1 + f * w); # g, dry weight at birth at f (remove d.V for wet weight)
    aT.b <- t.b/ k.M/ TC.ab;           # d, age at birth at f and T

    # puberty
    L.p <- L.m * l.p;                  # cm, structural length at puberty at f
    Lw.p <- L.p/ del.M;                # cm, physical length at puberty at f
    Wd.p <- L.p^3 * d.V * (1 + f * w); # g, dry weight at puberty (remove d.V for wet weight)
    aT.p <- t.p/ k.M/ TC.ap;           # d, age at puberty at f and T

    # ultimate
    l.i <- f - l.T;                    # -, scaled ultimate length at f
    L.i <- L.m * l.i;                  # cm, ultimate structural length at f
    Lw.i <- L.i/ del.M;                # cm, ultimate physical length at f
    Wd.i <- L.i^3 * d.V * (1 + f * w); # g, ultimate dry weight (remove d.V for wet weight)

    # reproduction
    pars.R <- c(kap = kap, kapR = kap.R, g = g, kJ = k.J, kM = k.M, LT = L.T, v = v, UHb = U.Hb, UHp = U.Hp); # compose parameter vector at T
    RT.i <- TC.Ri * reprod_rate(L.i, f, pars.R)[[1]];             # #/d, ultimate reproduction rate at T

    # life span
    pars.tm <- c(g = g, lT = l.T, ha = h.a/ k.M^2, sG = s.G);  # compose parameter vector at T.ref
    t.m <- get_tm_s(pars.tm, f, l.b)[[1]]; # -, scaled mean life span at T.ref
    aT.m <- t.m/ k.M/ TC.am;               # d, mean life span at T

    # pack to output
    # the names of the fields in the structure must be the same as the data names in the mydata file
    prdData <- list();
    prdData$ab <- aT.b;
    prdData$ap <- aT.p;
    prdData$am <- aT.m;
    prdData$Lb <- Lw.b;
    prdData$Lp <- Lw.p;
    prdData$Li <- Lw.i;
    prdData$Wdb <- Wd.b;
    prdData$Wdp <- Wd.p;
    prdData$Wdi <- Wd.i;
    prdData$Ri <- RT.i;

    # uni-variate data

    # time-length
    f <- f.tL; pars.lb <- c(g = g, k = k, vHb = v.Hb);          # compose parameters
    ir.B <- 3/ k.M + 3 * f * L.m/ v; r.B <- 1/ ir.B;            # d, 1/von Bert growth rate
    Lw.i <- (f * L.m - L.T)/ del.M;                             # cm, ultimate physical length at f
    Lw.b <- get_lb(pars.lb, f)[[1]] * L.m/ del.M;                    # cm, physical length at birth at f
    ELw <- Lw.i - (Lw.i - Lw.b) * exp( - TC.tL * r.B * tL[,1]); # cm, expected physical length at time
    #
    # length-weight
    EWw <- (LW[,1] * del.M)^3 * (1 + f * w);                    # g, expected wet weight at time

    # pack to output
    # the names of the fields in the structure must be the same as the data names in the mydata file
    prdData$tL <- ELw;
    prdData$LW <- EWw;

    return(list(prdData, info <- 1))
  }))))

}
