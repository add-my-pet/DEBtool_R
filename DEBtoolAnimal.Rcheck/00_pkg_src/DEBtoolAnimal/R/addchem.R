#' Sets chemical parameters and text for units and labels
#'
#' @description Sets chemical parameters and text for units and labels
#' @family add-my-pet auxiliary functions
#' @param par data frame with parameter values
#' @param units data frame with parameter units
#' @param label data frame with parameter labels
#' @param free data frame with information on which parameter to free or fix
#' @param phylum string with species phylum
#' @param class string with species class
#' @details Calls get_d_V to set specific density of structure. For a specific density of wet mass of 1 g/cm^3,
#' a specific density of d_E = d_V = 0.1 g/cm^3 means a dry-over-wet weight ratio of 0.1
#' @return list with updated par, units, label and free
#' @examples pars_init_my_pet(metaData)
#' @export
addchem <- function(par, units, label, free, phylum, class){

  # specific densities
  #   set specific densites using the pet's taxonomy
  list[d.V,] = get_d_V(phylum, class); # see comments on section 3.2.1 of DEB3
  par$d.X = d.V;     free$d.X = 0;  units$d.X = "g/cm^3"; label$d.X = "specific density of food";
  par$d.V = d.V;     free$d.V = 0;  units$d.V = "g/cm^3"; label$d.V = "specific density of structure";
  par$d.E = d.V;     free$d.E = 0;  units$d.E = "g/cm^3"; label$d.E = "specific density of reserve";
  par$d.P = d.V;     free$d.P = 0;  units$d.P = "g/cm^3"; label$d.P = "specific density of faeces";

  # chemical potentials from Kooy2010 Tab 4.2
  par$mu.X = 525000; free$mu.X = 0;  units$mu.X = "J/ mol"; label$mu.X = "chemical potential of food";
  par$mu.V = 500000; free$mu.V = 0;  units$mu.V = "J/ mol"; label$mu.V = "chemical potential of structure";
  par$mu.E = 550000; free$mu.E = 0;  units$mu.E = "J/ mol"; label$mu.E = "chemical potential of reserve";
  par$mu.P = 480000; free$mu.P = 0;  units$mu.P = "J/ mol"; label$mu.P = "chemical potential of faeces";

  # chemical potential of minerals
  par$mu.C = 0;      free$mu.C = 0;  units$mu.C = "J/ mol"; label$mu.C = "chemical potential of CO2";
  par$mu.H = 0;      free$mu.H = 0;  units$mu.H = "J/ mol"; label$mu.H = "chemical potential of H2O";
  par$mu.O = 0;      free$mu.O = 0;  units$mu.O = "J/ mol"; label$mu.O = "chemical potential of O2";
  par$mu.N = 0;      free$mu.N = 0;  units$mu.N = "J/ mol"; label$mu.N = "chemical potential of NH3";

  # chemical indices for water-free organics from Kooy2010 Fig 4.15 (excluding contributions of H and O atoms from water)
  par$n.CX = 1;      free$n.CX = 0;  units$n.CX = "-"; label$n.CX = "chem. index of carbon in food"; # C/C = 1 by definition
  par$n.HX = 1.8;    free$n.HX = 0;  units$n.HX = "-"; label$n.HX = "chem. index of hydrogen in food";
  par$n.OX = 0.5;    free$n.OX = 0;  units$n.OX = "-"; label$n.OX = "chem. index of oxygen in food";
  par$n.NX = 0.15;   free$n.NX = 0;  units$n.NX = "-"; label$n.NX = "chem. index of nitrogen in food";
  #
  par$n.CV = 1;      free$n.CV = 0;  units$n.CV = "-"; label$n.CV = "chem. index of carbon in structure"; # n.CV = 1 by definition
  par$n.HV = 1.8;    free$n.HV = 0;  units$n.HV = "-"; label$n.HV = "chem. index of hydrogen in structure";
  par$n.OV = 0.5;    free$n.OV = 0;  units$n.OV = "-"; label$n.OV = "chem. index of oxygen in structure";
  par$n.NV = 0.15;   free$n.NV = 0;  units$n.NV = "-"; label$n.NV = "chem. index of nitrogen in structure";
  #
  par$n.CE = 1;      free$n.CE = 0;  units$n.CE = "-"; label$n.CE = "chem. index of carbon in reserve";   # n.CE = 1 by definition
  par$n.HE = 1.8;    free$n.HE = 0;  units$n.HE = "-"; label$n.HE = "chem. index of hydrogen in reserve";
  par$n.OE = 0.5;    free$n.OE = 0;  units$n.OE = "-"; label$n.OE = "chem. index of oxygen in reserve";
  par$n.NE = 0.15;   free$n.NE = 0;  units$n.NE = "-"; label$n.NE = "chem. index of nitrogen in reserve";
  #
  par$n.CP = 1;      free$n.CP = 0;  units$n.CP = "-"; label$n.CP = "chem. index of carbon in faeces";    # n.CP = 1 by definition
  par$n.HP = 1.8;    free$n.HP = 0;  units$n.HP = "-"; label$n.HP = "chem. index of hydrogen in faeces";
  par$n.OP = 0.5;    free$n.OP = 0;  units$n.OP = "-"; label$n.OP = "chem. index of oxygen in faeces";
  par$n.NP = 0.15;   free$n.NP = 0;  units$n.NP = "-"; label$n.NP = "chem. index of nitrogen in faeces";

  # chemical indices for minerals from Kooy2010
  par$n.CC = 1;   free$n.CC = 0;  units$n.CC = "-"; label$n.CC = "chem. index of carbon in CO2";
  par$n.HC = 0;   free$n.HC = 0;  units$n.HC = "-"; label$n.HC = "chem. index of hydrogen in CO2";
  par$n.OC = 2;   free$n.OC = 0;  units$n.OC = "-"; label$n.OC = "chem. index of oxygen in CO2";
  par$n.NC = 0;   free$n.NC = 0;  units$n.NC = "-"; label$n.NC = "chem. index of nitrogen in CO2";
  #
  par$n.CH = 0;   free$n.CH = 0;  units$n.CH = "-"; label$n.CH = "chem. index of carbon in H2O";
  par$n.HH = 2;   free$n.HH = 0;  units$n.HH = "-"; label$n.HH = "chem. index of hydrogen in H2O";
  par$n.OH = 1;   free$n.OH = 0;  units$n.OH = "-"; label$n.OH = "chem. index of oxygen in H2O";
  par$n.NH = 0;   free$n.NH = 0;  units$n.NH = "-"; label$n.NH = "chem. index of nitrogen in H2O";
  #
  par$n.CO = 0;   free$n.CO = 0;  units$n.CO = "-"; label$n.CO = "chem. index of carbon in O2";
  par$n.HO = 0;   free$n.HO = 0;  units$n.HO = "-"; label$n.HO = "chem. index of hydrogen in O2";
  par$n.OO = 2;   free$n.OO = 0;  units$n.OO = "-"; label$n.OO = "chem. index of oxygen in O2";
  par$n.NO = 0;   free$n.NO = 0;  units$n.NO = "-"; label$n.NO = "chem. index of nitrogen in O2";
  #
  par$n.CN = 0;   free$n.CN = 0;  units$n.CN = "-"; label$n.CN = "chem. index of carbon in NH3";
  par$n.HN = 3;   free$n.HN = 0;  units$n.HN = "-"; label$n.HN = "chem. index of hydrogen in NH3";
  par$n.ON = 0;   free$n.ON = 0;  units$n.ON = "-"; label$n.ON = "chem. index of oxygen in NH3";
  par$n.NN = 1;   free$n.NN = 0;  units$n.NN = "-"; label$n.NN = "chem. index of nitrogen in NH3";

  return(list(par = par, units = units, label = label, free = free))
}
