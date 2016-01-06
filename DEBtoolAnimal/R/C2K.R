#' Conversion of Celsius to Kelvin
#'
#' @description Converts temperature in degrees Celsius to Kelvin
#' @family miscellaneous functions
#' @param C numeric temperature in degrees Celsius
#' @return temperature in Kelvin
#' @examples
#' C2K(20)
#' @export
C2K <- function(C) {
  return(C + 273.15)   # K, temperature in Kelvin
}
