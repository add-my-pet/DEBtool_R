#' Checks if parameters allow for reaching birth in the standard DEB model
#'
#' @description Checks if parameters allow for reaching birth in the standard DEB model
#' @family filter functions
#' @param g energy investment ratio
#' @param k ratio of maturity and somatic maintenance rate coeff
#' @param vHb scaled maturity volume at birth
#' @param f functional response (default 1)
#' @return info, indicator equals 1 if reaches birth, 0 otherwise
#' @examples reach_birth(g = 10, k = 1, vHb = 0.5)
#' @export
reach_birth <- function(g, k, vHb, f = 1){

  list[l.b, info] <- get_lb(c(g = g, k = k, vHb = vHb), f)

  if(l.b >= 1 || !info)
    return(info <- 0)

  if(k * vHb >= f/ (g + f) * l.b^2 * (g + l.b))
    return(info <- 0)

  return(info <- 1)

}
