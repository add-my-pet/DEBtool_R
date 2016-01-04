#' Adds pseudodata information into inputed data structures
#'
#' @description Adds the pseudodata information and weights for purposes of the regression
#' @family add-my-pet auxiliary functions
#' @param data structure with data
#' @param units structure with data units
#' @param label structure with data labels
#' @param weights structure with weights
#' @return structures with data, units, labael and weights
#' @examples list[data, units, label, weight] <- addpseudodata();
#' @export
addpseudodata <- function(data = list(), units = list(), label = list(), weights = list()){
  # set pseudodata
  data$psd$v <- 0.02;     units$psd$v <- "cm/d";       label$psd$v <- "energy conductance";
  data$psd$kap <- 0.8;    units$psd$kap <- "-";        label$psd$kap <- "allocation fraction to soma";
  data$psd$kap.R <- 0.95; units$psd$kap.R <- "-";      label$psd$kap.R <- "reproduction efficiency";
  data$psd$p.M <- 18;     units$psd$p.M <- "J/d.cm^3"; label$psd$p.M <- "vol-spec som maint";
  data$psd$p.T <- 0;      units$psd$p.T <- "J/d.cm^2"; label$psd$p.T <- "surf-spec som maint";
  data$psd$k.J <- 0.002;  units$psd$k.J <- "1/d";      label$psd$k.J <- "maturity maint rate coefficient";
  data$psd$kap.G <- 0.8;  units$psd$kap.G <- "-";      label$psd$kap.G <- "growth efficiency";

  # set weights
  weights$psd <- setweights(data$psd);
  weights$psd$v     <- 0.1 * weights$psd$v;
  weights$psd$kap   <- 0.1 * weights$psd$kap;
  weights$psd$kap.R <- 0.1 * weights$psd$kap.R;
  weights$psd$p.M   <- 0.1 * weights$psd$p.M;
  weights$psd$p.T   <- 0.1 * weights$psd$p.T;
  weights$psd$k.J   <- 0.1 * weights$psd$k.J;
  weights$psd$kap.G <- 0.1 * weights$psd$kap.G;

  weights$psd$kap.G <- 200 * weights$psd$kap.G;   # more weight to kap.G

  return(list(data, units, label, weights))

}
