#' Calculates least squares estimates using Nelder Mead's simplex method using a filter
#'
#' @description Calculates least squares estimates using Nelder Mead's simplex method using a filter
#' @family add-my-pet functions
#' @param func character string with name of user-defined function
#' @param par list with parameters
#' @param data list with data
#' @param auxData list with auxiliary data
#' @param weights list with weights
#' @param filternm character string with name of user-defined filter function
#' @return list with list with parameters resulting from estimation procedure (par)
#' and indicator 1 if convergence has been successful or 0 otherwise (info)
#' @export
petregr_f <- function(func, par, data, auxData, weights, filternm){

  # option settings
  info <- 1; # initiate info setting

  # prepare variable
  #   st: structure with dependent data values only
  st <- data;
  nst <- fieldnmnst(st); # nst: number of data sets

  for(currentData in nst) {   # makes st only with dependent variables
    auxVar <- st[[currentData]];   # data in field currentData
    k <- size(auxVar, 2);
    if(k > 1)
      st <- auxvar[,2];
  }

  # Y: vector with all dependent data
  # W: vector with all weights
  Y <- struct2vector(st, nm);
  W <- struct2vector(weights, nm);




  return(list(par = q, info))
}



