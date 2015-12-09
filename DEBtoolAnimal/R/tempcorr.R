#' Temperature correction
#'
#' @description Calculates the factor with which physiological rates should be multiplied to go from a reference temperature to a given temperature
#' @family miscelaneous functions
#' @param Temp vector with temperatures (in Kelvin)
#' @param T_1 scalar with reference temperature (in Kelvin)
#' @param T_A scalar with Arrhenius temperature (in Kelvin)
#' @param T_L optional scalar with lower boundary of temperature range (in Kelvin)
#' @param T_AL optional scalar with Arrhenius temperature for lower boundary of temperature range (in Kelvin)
#' @param T_H optional scalar with upper boundary of temperature range (in Kelvin)
#' @param T_AH optional scalar with Arrhenius temperature for upper boundary of temperature range (in Kelvin)
#' @return vector with temperature correction factors that affect all rates
#' @details Temperature impacts metabolic rates. This impact, in its most simplest way (1 parameter), is modeled by multiplying all the time-dependent parameters by a correction factor:
#' \deqn{\exp\left(\frac{T_A}{T_1} - \frac{T_A}{T}\right)}{exp(T_A/ T_1 - T_A/ T)}
#' For a more detailed modeling one can multiply with an extra fraction \eqn{s(T_1)/s(T)} with (3 parameters):
#' \deqn{s(T) = 1 + \exp\left(\frac{T_{AL}}{T} - \frac{T_{AL}}{T_L}\right}{s(T) = 1 + exp(T_AL/ T - T_AL/ T_L))}
#' or (5 parameters)
#' \deqn{s(T) = 1 + \exp\left(\frac{T_{AL}}{T} - \frac{T_{AL}}{T_L}\right + \exp\left(\frac{T_{AH}}{T_H} - \frac{T_{AH}}{T}\right}{s(T) = 1 + exp(T_AL/ T - T_AL/ T_L) + exp(T_AH/ T_H - T_AH/ T)}
#' @examples tempcorr(c(330, 331, 332), 320, T_A = 12000, T_L = 277, T_H = 331, T_AL = 20000, T_AH = 190000)
#' @export
tempcorr <- function(Temp, T_1, T_A, T_L = NA, T_AL = NA, T_H = NA, T_AH = NA){
  termT  <- 0
  termT1 <- 0

  if (!is.na(T_L * T_AL)){
    termT <- termT + exp(T_AL/ Temp - T_AL/ T_L)
    termT1 <- termT1 + exp(T_AL/ T_1 - T_AL/ T_L)
  }

  if (!is.na(T_H * T_AH)){
    termT <- termT + exp(T_AH/ T_H - T_AH/ Temp)
    termT1 <- termT1 + exp(T_AH/ T_H - T_AH/ T_1)
  }

  exp(T_A/ T_1 - T_A/ Temp) * (1 + termT1)/  (1 + termT)
}
