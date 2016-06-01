#' Predicts pseudodata values
#'
#' @description Adds pseudodata predictions into predictions structure
#' @family add-my-pet auxiliary functions
#' @param par data frame with parameter values
#' @param data data frame with data values
#' @param prdData data frame with prediction values
#' @return structure with pseudodata predictions
#' @export
predict_pseudodata <- function(par, data, prdData){

  nm <- fieldnm_wtxt(data, "psd");

  if(length(nm) > 0)
    # unpack coefficients
    cPar <- parscomp(par);  allPar <- par;
    fieldNames <- names(cPar);
    for(currentCmpPar in fieldNames)
      allPar[[currentCmpPar]] <- cPar[[currentCmpPar]];

  varnm <- names(data[[nm[[1]]]]);

  # adds pseudodata predictions to structure
  prdData$psd <- list();
  for(currentPsd in varnm)
    prdData$psd[[currentPsd]] <- allPar[[currentPsd]];

  return(prdData)
}
