#' Prints an explanation of the filter flag onto the screen
#'
#' @description Prints an explanation to the screen according to the flag produced by a filter.
#' Meant to be run in the estimation procedure for the seed parameter set
#' @family add-my-pet auxiliary functions
#' @param flag integer with code from filter
#' @examples print_filterflag(3)
#' @export
print_filterflag <- function(flag){

  cat(
    switch(toString(flag),
         "1" = "One or more parameters are negative \n",
         "2" = "kappa, f or one of the fractions is bigger than 1 \n",
         "3" = "growth efficiency is bigger than 1 \n",
         "4" = "maturity levels are not in the correct order \n",
         "5" = "puberty or birth cannot be reached \n"
    )
  )

}
