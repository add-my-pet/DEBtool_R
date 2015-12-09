#' Conversion of Kelvin to Celsius
#'
#' @description Converts temperature in Kelvin to degrees Celsius
#' @family miscellaneous functions
#' @param K numeric temperature in degrees Kelvin
#' @return temperature in Kelvin
#' @examples
#' K2C(293.15)
#' @export
K2C <- function(K) {
  return(K - 273.15)   # C, temperature in Celsius
}
