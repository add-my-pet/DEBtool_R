#' Prints data of a species to screen
#'
#' @description Prints data of a species to screen
#' @family add-my-pet functions
#' @param data list with data values
#' @param txtData list with text info on data
#' @param prdData list with prediction values
#' @param RE list with relative errors
#' @export
printprd <- function(data, txtData, prdData, RE){

  nm <- fieldnmnst(data);
  dtsets <- names(data);

  cat("data and predictions (relative error) \n");
  j <- 1;
  for(currentData in nm) {
    tempData <- data[[currentData]];
    k <- dim(as.matrix(tempData))[2]; # number of data points per set
    if(k == 1)
      if(!is.list(currentData)) {
        tempUnit <- txtData$units[[currentData]];
        tempLabel <- txtData$label[[currentData]];
        tempPrdData <- prdData[[currentData]];
        cat(tempData, " ", tempPrdData, " (", RE[j], ") ", currentData, ", ", tempUnit, ", ", tempLabel, "\n", sep = "");
      } else {
        cat(tempData, " ", tempPrdData, " (", RE[j], ") ", currentData, ", ", tempUnit, ", ", tempLabel, "\n", sep = "");
      }
      else {
        cat("see figure (", RE[j], ") ", currentData, ", ", txtData$label[[currentData]][[2]], " vs. ", txtData$label[[currentData]][[1]], "\n", sep = "");
      }
    j <- j + 1;

  }
}
