#' Kelvin form Celsius
#' 
#' converts degrees Celsius to Kelvin
#' 
#' @description Obtains Kelvin from temperatures defined in Celsius
#' @param C: numeric temperature in degrees Celsius
#' @return K: temperature in Kelvin
#' @example C2K(20) 
#' @export
C2K <- function(C) {
  K <- C + 273.15   # K, temperature in Kelvin
  return(K)
}
