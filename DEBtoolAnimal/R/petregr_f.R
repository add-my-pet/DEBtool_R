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




  return(list(par = q, info))
}


cov_rules_1species <- function(par, i) par


