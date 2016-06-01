#' Concatenates predict files for several species
#'
#' @description Concatenates predict files for several species
#' @family add-my-pet functions
#' @param parGrp data frame with parameter values of the group
#' @param data data frame with data values
#' @param auxData data frame with auxiliary data values
#' @return structure with prdData and prdInfo for several pets
#' @export
predict_pets <- function(parGrp, data, auxData){

  info <- 0;

  cov.rulesnm <- paste("cov_rules_", cov_rules, sep = "")

  # produce pars for species and predict
  prdData <- list();
  i <- 1;
  for(currentPet in pets) {
    # for the case with no zoom factor transformation
    par <- do.call(cov.rulesnm, list(parGrp, i));
    list[prdData[[currentPet]], info] <- do.call(paste("predict_", currentPet, sep = ""), list(par, data[[currentPet]], auxData[[currentPet]]));
    if(!info)
      return(prdData = NA, info);

    # predict pseudodata
    if(pseudodata_pets == 0) # option of estim
      prdData[[currentPet]] <- predict_pseudodata(par, data[[currentPet]], prdData[[currentPet]]);

    i <- i + 1;
  }

  if(pseudodata_pets == 1)
    # predicts pseudodata
    prdData <- predict_pseudodata(do.call(covRulesnm, list(parGrp, 1)), data, prdData);

  info <- 1;

  return(list(prdData = prdData, info = info))
}

cov_rules_1species <- function(par, i) par
