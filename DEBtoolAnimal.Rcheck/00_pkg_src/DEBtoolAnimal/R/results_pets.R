#' Prints results of estimation
#'
#' @description Prints the results of the esimation procedure in the screen, .mat file and makes figures of graphs
#' @family add-my-pet functions
#' @param par data frame with parameter values
#' @param metaPar data frame with metainformation on models
#' @param txtPar data frame with information on parameters
#' @param data data frame with data values
#' @param auxData data frame with auxiliary data values
#' @param metaData data frame with metainformation on the entry
#' @param txtData data frame with infromation on data
#' @param weights data frame with values of weights
#' @examples results_pets(par, metaPar, txtPar, data, auxData, metaData, txtData, weights)
#' @export
results_pets <- function(par, metaPar, txtPar, data, auxData, metaData, txtData, weights){

  petsNumber <- length(pets)

  cov.rulesnm <- paste("cov_rules_", metaPar$covRules, sep = "")

  weightsMRE <- weights;  # define a weights structure with weight 1 for every data point and 0 for the pseudodata
  for(currentPet in pets) {   # makes st only with dependent variables
    if("psd" %in% weights[[currentPet]]) {
      psdSets = names(weights[[currentPet]]$psd);
      for(currentPsd in psdSets)
        weightsMRE[[currentPet]]$psd[[currentPsd]] = 0 * weights[[currentPet]]$psd[[currentPsd]];
    }
  }




}

cov_rules_1species <- function(par, i) par
