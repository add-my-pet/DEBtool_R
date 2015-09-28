## C2K
# computes Kelvin from Celsius

C2K <- function(C) {
  ## Description
  # Obtains Kelvin from temperatures defined in Celsius
  #
  # Input
  #
  # * C: scalar or matrix in temperatures in degrees Celsius
  #  
  # Output
  #
  # * K: temperature(s) in Kelvin
  
  ## Example 
  # C2K(20)
  
  
  K <- C + 273.15   # K, temperature in Kelvin
  return(K)
}
