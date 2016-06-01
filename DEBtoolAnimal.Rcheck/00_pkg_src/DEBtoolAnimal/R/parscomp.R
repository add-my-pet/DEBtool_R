#' Computes compound parameters from primary parameters
#'
#' @description Computes compound parameters from primary parameters that are frequently used
#' @family add-my-pet auxiliary functions
#' @param par data frame with parameter values
#' @return list with compound parameters
#' @examples parscomp(par)
#' @export
parscomp <- function(par){

  with(as.list(par), {

    if(!exists("p.Am"))
      p.Am = z * p.M/ kap         # J/d.cm^2, {p_Am} spec assimilation flux

    ##              X     V     E     P
    n.O <- cbind(c(n.CX, n.CV, n.CE, n.CP),  # C/C, equals 1 by definition
                 c(n.HX, n.HV, n.HE, n.HP),  # H/C  these values show that we consider dry-mass
                 c(n.OX, n.OV, n.OE, n.OP),  # O/C
                 c(n.NX, n.NV, n.NE, n.NP))  # N/C

    ##              C     H     O     N
    n.M <- cbind(c(n.CC, n.CH, n.CO, n.CN),  # CO2
                 c(n.HC, n.HH, n.HO, n.HN),  # H2O
                 c(n.OC, n.OH, n.OO, n.ON),  # O2
                 c(n.NC, n.NH, n.NO, n.NN))  # NH3

    # -------------------------------------------------------------------------
    # Molecular weights:
    w.O <- n.O %*% c(12, 1, 16, 14) # g/mol, mol-weights for (unhydrated)  org. compounds
    w.X <- w.O[1];
    w.V <- w.O[2];
    w.E <- w.O[3];
    w.P <- w.O[4];

    # -------------------------------------------------------------------------
    # Conversions and compound parameters cPar
    M.V     <- d.V/ w.V;            # mol/cm^3, [M_V], volume-specific mass of structure
    y.V.E   <- mu.E * M.V/ E.G;     # mol/mol, yield of structure on reserve
    y.E.V   <- 1/ y.V.E;            # mol/mol, yield of reserve on structure
    k.M     <- p.M/ E.G;            # 1/d, somatic maintenance rate coefficient
    k       <- k.J/ k.M;            # -, maintenance ratio
    E.m     <- p.Am/ v;             # J/cm^3, [E_m], reserve capacity
    m.Em    <- y.E.V * E.m / E.G;   # mol/mol, reserve capacity
    g       <- E.G/ kap/ E.m ;      # -, energy investment ratio
    L.m     <- v/ k.M/ g;           # cm, maximum length
    L.T     <- p.T/ p.M ;           # cm, heating length (also applies to osmotic work)
    l.T     <- L.T/ L.m;            # - , scaled heating length
    w       <- m.Em * w.E * d.V/ d.E/ w.V; # -, \omega, contribution of ash free dry mass of reserve to total ash free dry biomass
    J.E.Am  <- p.Am/ mu.E;          # mol/d.cm^2, {J_EAm}, max surface-spec assimilation flux

    if(exists("kap.X")) {
      y.E.X  <- kap.X * mu.X/ mu.E; # mol/mol, yield of reserve on food
      y.X.E  <- 1/ y.E.X;           # mol/mol, yield of food on reserve
      p.Xm   <- p.Am/ kap.X;        # J/d.cm^2, max spec feeding power
      J.X.Am <- y.X.E * J.E.Am;     # mol/d.cm^2, {J_XAm}, max surface-spec feeding flux
    }

    if(exists("kap.P")) {
      y.P.X  <- kap.P * mu.X/ mu.P; # mol/mol, yield of faeces on food
      y.X.P  <- 1/ y.P.X;           # mol/mol, yield of food on faeces
    }

    if(exists("kap.X") && exists("kap.P")) {
      y.P.E  <- y.P.X/ y.E.X;       # mol/mol, yield of faeces on reserve
      #  Mass-power couplers
      eta.XA <- y.X.E/ mu.E;        # mol/J, food-assim energy coupler
      eta.PA <- y.P.E/ mu.E;        # mol/J, faeces-assim energy coupler
      eta.VG <- y.V.E/ mu.E;        # mol/J, struct-growth energy coupler
      eta.O  <- cbind(c(-eta.XA,       0,       0),    # mol/J, mass-energy coupler
                      c(      0,       0,  eta.VG),    # used in: J_O = eta_O * p
                      c( 1/mu.E, -1/mu.E, -1/mu.E),
                      c( eta.PA,       0,       0));
    }

    J.E.M   <- p.M/ mu.E;          # mol/d.cm^3, [J_EM], volume-spec somatic  maint costs
    J.E.T   <- p.T/ mu.E;          # mol/d.cm^2, {J_ET}, surface-spec somatic  maint costs
    j.E.M   <- k.M * y.E.V;        # mol/d.mol, mass-spec somatic  maint costs
    j.E.J   <- k.J * y.E.V;        # mol/d.mol, mass-spec maturity maint costs
    kap.G   <- mu.V * M.V/ E.G;    # -, growth efficiency
    E.V     <- d.V * mu.V/ w.V;    # J/cm^3, [E_V] volume-specific energy of structure

    if(exists("F.m"))
      K <- J.X.Am/ F.m;        # c-mol X/l, half-saturation coefficient


    # -------------------------------------------------------------------------
    # pack output:
    cPar <- list(p.Am = p.Am, n.O = n.O, n.M = n.M, w.X = w.X, w.V = w.V, w.E = w.E, w.P = w.P, M.V = M.V,
                 y.V.E = y.V.E, y.E.V = y.E.V, k.M = k.M, k = k, E.m = E.m, m.Em = m.Em, g = g, L.m = L.m,
                 L.T = L.T, l.T = l.T, w = w, J.E.Am = J.E.Am, J.E.M = J.E.M, J.E.T = J.E.T, j.E.M = j.E.M,
                 j.E.J = j.E.J, kap.G = kap.G, E.V = E.V)

    # -------------------------------------------------------------------------
    # add the Scaled maturity maturity levels:
    parNames <- names(par);
    matInd <- gsub("E.H", "", parNames[grep("E.H", parNames)]);  # maturity levels' indices

    for(i in matInd) {
      cPar[paste("M.H", i, sep = "")] <- par[[paste("E.H", i, sep = "")]]/ mu.E;               # mmol, maturity at level i
      cPar[paste("U.H", i, sep = "")] <- par[[paste("E.H", i, sep = "")]]/ p.Am;               # cm^2 d, scaled maturity at level i
      cPar[paste("V.H", i, sep = "")] <- cPar[[paste("U.H", i, sep = "")]]/ (1 - kap);         # cm^2 d, scaled maturity at level i
      cPar[paste("v.H", i, sep = "")] <- cPar[[paste("V.H", i, sep = "")]] * g^2 * k.M^3/ v^2; # -, scaled maturity density at level i
      cPar[paste("u.H", i, sep = "")] <- cPar[[paste("U.H", i, sep = "")]] * g^2 * k.M^3/ v^2; # -, scaled maturity density at level i
    }

    if(exists("kap.X")) {
      cPar$y.E.X  <- y.E.X
      cPar$y.X.E  <- y.X.E
      cPar$p.Xm   <- p.Xm
      cPar$J.X.Am <- J.X.Am
    }

    if(exists("kap.P")) {
      cPar$y.P.X  <- y.P.X
      cPar$y.X.P  <- y.X.P
    }

    if(exists("kap.X") && exists("kap.P")) {
      cPar$y.P.E  <- y.P.E
      cPar$eta.XA <- eta.XA
      cPar$eta.PA <- eta.PA
      cPar$eta.VG <- eta.VG
      cPar$eta.O  <- eta.O
    }

    if(exists("F.m"))
      cPar$K <- K

    return(cPar)
  })
}



