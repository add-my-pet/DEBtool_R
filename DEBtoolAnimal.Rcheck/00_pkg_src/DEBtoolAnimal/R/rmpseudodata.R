#' Removes pseudodata information from inputed data structures
#'
#' @description Removes pseudodata information from inputed data structures
#' @family add-my-pet auxiliary functions
#' @param data structure with "psd" field to be removed
#' @return structure with "psd" field removed
#' @examples data <- rmpseudodata(data)
#' @export
rmpseudodata <- function(data = list()){

  nmpsd <- fieldnm_wtxt(data, "psd");

  if(length(nmpsd) > 0) {
    for(currentPsd in nmpsd)
      eval(parse(text = paste("data$", currentPsd, " <- NULL", sep = "")))
  }

  return(data)
}
