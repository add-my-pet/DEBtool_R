#' Particular incomplete beta function
#'
#' @description particular incomplete beta function:
#   B_x1(4/3,0) - B_x0(4/3,0) = \int_x0^x1 t^(4/3-1) (1-t)^(-1) dt.
#' @family miscellaneous functions
#' @param x0 scalar with lower boundary for integration
#' @param x1 scalar with upper boundary for integration
#' @return scalar with particular incomple beta function
#' @details Computes
#' \deqn{B_{x_1}\left(\frac{4}{3},0\right) - B_{x_0}\left(\frac{4}{3},0\right) = \int_{x_0}^{x_1} t^{4/3-1} (1-t)^{-1} dt}{B_x1(4/3,0) - B_x0(4/3,0) = integral (from x0 to x1) t^(4/3-1) (1-t)^(-1) dt}
#' To be used in the computation of the age at birth (or related quantities) for an egg.
#' @examples
#' beta0(0.1, 0.2)
#' @export
beta0 <- function(x0, x1){

  if (x0 < 0 || x0 >= 1 || x1 < 0 || x1 >= 1){
    print('Warning from beta0: argument values outside (0,1) ')
    return(NaN)
    #stop("argument values outside (0,1)")
  }

  n0 <- length(x0)
  n1 <- length(x1)
  if (n0 != n1 && n0 != 1 && n1 != 1) {
    stop("argument sizes do not match")
  }

  x03 <- x0 ^ (1/ 3)
  x13 <- x1 ^ (1/ 3)
  a3 <- sqrt(3)
  3 * (x03 - x13) + a3 * atan((1 + 2 * x13)/ a3) - a3 * atan((1 + 2 * x03)/ a3) +
    log((x03 - 1)/ (x13 - 1)) + log((1 + x13 + x13 ^ 2)/ (1 + x03 + x03 ^ 2))/ 2
}
