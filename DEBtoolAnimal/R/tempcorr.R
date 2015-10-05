#' Conversion of Kelvin to Celsius
#'
#' @description Calculates the factor with which physiological rates should be multiplied to go from a reference temperature to a given temperature
#' @family miscelaneous functions
#' @param T vector with new temperatures
#' @param T_1 scalar with reference temperature
#' @param Tpars 1-, 3- or 5-vector with temperature parameters
#' @return vector with temperature correction factors that affect all rates
#' @examples tempcorr(c(330, 331, 332), 320, c(12000, 277, 318, 20000, 190000))
#' @export
tempcorr <- function(Temp, T_1, Tpars){
  nmPars <- length(Tpars)
  T_A <- Tpars[1]; # Arrhenius temperature
  switch(as.character(nmPars),
         "1" = exp(T_A/ T_1 - T_A/ Temp),
         "3" = {T_L  <- Tpars[2]  # Lower temp boundary
                T_AL <- Tpars[3]  # Arrh. temp for lower boundary
                exp(T_A/ T_1 - T_A/ Temp) * (1 + exp(T_AL/ T_1 - T_AL/ T_L)) / (1 + exp(T_AL/ Temp - T_AL/ T_L)) },
         "5" = {T_L  <- Tpars[2]  # Lower temp boundary
                T_H  <- Tpars[3]  # Upper temp boundary
                T_AL <- Tpars[4]  # Arrh. temp for lower boundary
                T_AH <- Tpars[5]  # Arrh. temp for upper boundary
                exp(T_A/ T_1 - T_A/ Temp) *
                  (1 + exp(T_AL/ T_1 - T_AL/ T_L) + exp(T_AH/ T_H - T_AH/ T_1))/
                  (1 + exp(T_AL/ Temp - T_AL/ T_L) + exp(T_AH/ T_H - T_AH/ Temp)) },
         stop("The number of parameters for the temperature correction must be 1, 3 or 5.")
         )
}
