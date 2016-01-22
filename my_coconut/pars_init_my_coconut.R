#' Sets (initial values for) parameters
#'
#' @description Sets (initial values for) parameters$ Meant to be a template in add-my-pet
#' @family add-my-pet template functions
#' @param metaData data frame with info about this entry (needed for names of phylum and class to get d_V)
#' @return list with par (with values of parameters), metaPar (with information on metaparameters) and txtPar (with information on parameters)
#' @examples pars_init_my_pet(metaData)
#' @export
pars_init_coconut <- function(metaData){

  par <- list(); free <- list(); units <- list(); label <- list();
  metaPar <- list(); txtPar <- list();

  # parameters: initial values at reference temperature T_ref
  # model type: see online manual for explanation and alternatives
  # be aware that each model type is associated with a specific list of core
  # primary parameters$ Those listed here are for model "std". See the manual
  # for parameters associated with the other model types.
  metaPar$model <- "std";

  # edit the values below such that predictions are not too far off;
  # the values must be set in the standard DEB units:
  #   d for time; cm for length; J for energy; K for temperature

  # reference parameter (not to be changed)
  par$T.ref <- C2K(20); free$T.ref <- 0; units$T.ref <- "K";        label$T.ref <- "Reference temperature";

  ## core primary parameters
  par$z     <- 1;     free$z     <- 1;   units$z     <- "-";        label$z     <- "zoom factor"; #for z = 1: L_m = 1 cm
  par$F.m   <- 6.5;   free$F.m   <- 0;   units$F.m   <- "l/d.cm^2"; label$F.m   <- "{F_m}, max spec searching rate";
  par$kap.X <- 0.8;   free$kap.X <- 0;   units$kap.X <- "-";        label$kap.X <- "digestion efficiency of food to reserve";
  par$kap.P <- 0.1;   free$kap.P <- 0;   units$kap.P <- "-";        label$kap.P <- "faecation efficiency of food to faeces";
  par$v     <- 0.02;  free$v     <- 1;   units$v     <- "cm/d";     label$v     <- "energy conductance";
  par$kap   <- 0.8;   free$kap   <- 1;   units$kap   <- "-";        label$kap   <- "allocation fraction to soma";
  par$kap.R <- 0.95;  free$kap.R <- 0;   units$kap.R <- "-";        label$kap.R <- "reproduction efficiency";
  par$p.M   <- 18;    free$p.M   <- 1;   units$p.M   <- "J/d.cm^3"; label$p.M   <- "[p_M], vol-spec somatic maint";
  par$p.T   <-  0;    free$p.T   <- 0;   units$p.T   <- "J/d.cm^2"; label$p.T   <- "{p_T}, surf-spec somatic maint";
  par$k.J   <- 0.002; free$k.J   <- 0;   units$k.J   <- "1/d";      label$k.J   <- "maturity maint rate coefficient";
  par$E.G   <- 2800;  free$E.G   <- 1;   units$E.G   <- "J/cm^3";   label$E.G   <- "[E_G], spec cost for structure";
  par$E.Hb  <- .275;  free$E.Hb  <- 1;   units$E.Hb  <- "J";        label$E.Hb  <- "maturity at birth";
  par$E.Hp  <- 50;    free$E.Hp  <- 1;   units$E.Hp  <- "J";        label$E.Hp  <- "maturity at puberty";
  par$h.a   <- 1e-6;  free$h.a   <- 1;   units$h.a   <- "1/d^2";    label$h.a   <- "Weibull aging acceleration";
  par$s.G   <- 1e-4;  free$s.G   <- 0;   units$s.G   <- "-";        label$s.G   <- "Gompertz stress coefficient";

  ## auxiliary parameters
  par$T.A   <- 8000;   free$T.A   <- 0;    units$T.A   <- "K";      label$T.A   <- "Arrhenius temperature";
  par$del.M <- 0.16;   free$del.M <- 1;    units$del.M <- "-";      label$del.M <- "shape coefficient";

  ## environmental parameters (temperatures are in data)
  par$f    <- 1.0;     free$f     <- 0;    units$f    <- "-";       label$f     <- "scaled functional response for 0-var data";
  par$f.tL <- 0.8;     free$f.tL  <- 1;    units$f.tL <- "-";       label$f.tL  <- "scaled functional response for 1-var data";

  ## set chemical parameters from Kooy2010
  #  don't change these values, unless you have a good reason
  list[par, units, label, free] = addchem(par, units, label, free, phylum = metaData$phylum, class = metaData$class);

  # if you do have a good reason you may overwrite any of the values here, but please provide an explanations.
  # For example:
  # par$d.V = 0.33;    # g/cm^3, specific density of structure, see ref bibkey
  # par$d.E = par$d.V; # g/cm^3, specific density of reserve
  #   or alternatively you might want to add a product D like:
  # par$d.D = 0.1;     free$d.D = 0;  units$d.D  = "g/cm^3"; label$d.D  = "specific density of product";
  # par$mu.D = 273210; free$mu.D = 0; units$mu.D = "J/mol";  label$mu.D = "chemical potential of product";
  # par$n.CD = 1;      free$n.CD = 0; units$n.CD = "-";      label$n.CD = "chem. index of carbon in product";
  # par$n.HD = 1.2;    free$n.HD = 0; units$n.HD = "-";      label$n.HD = "chem. index of hydrogen in product";
  # par$n.OD = 0.55;   free$n.OD = 0; units$n.OD = "-";      label$n.OD = "chem. index of oxygen in product";
  # par$n.ND = 0.1;    free$n.ND = 0; units$n.ND = "-";      label$n.ND = "chem. index of nitrogen in product";

  ## estimating chemical parameters (remove these remarks after editing the file)
  # in some cases the data allow for estimating some of the chemical
  # parameters. For example:
  # free$mu_V = 1;
  # free$d_E  = 1;

  ## Pack output:
  txtPar$units <- units; txtPar$label <- label; par$free <- free;

  return(list(par = par, metaPar = metaPar, txtPar = txtPar))
}
