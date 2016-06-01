#' Prints parameters of a species to screen
#'
#' @description Prints parameters of a species to screen
#' @family add-my-pet functions
#' @param par list with parameter values
#' @param txtPar list with text info on parameters
#' @export
printpar <- function(par, txtPar){

  cat("\nParameters \n");
  free <- par$free;                 # copy substructure
  parpl <- par; parpl$free <- NULL; # remove substructure free from par
  nm <- fieldnmnst(parpl);          # get number of parameter fields
  for(currentPar in nm) # scan parameter fields
    cat(currentPar, ", ", txtPar$units[[currentPar]], " ", parpl[[currentPar]], " ", free[[currentPar]], "\n", sep = "");
}
