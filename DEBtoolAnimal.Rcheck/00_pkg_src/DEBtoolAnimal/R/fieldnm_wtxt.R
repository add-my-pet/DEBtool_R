#' Searches for fields with a given name in a multilevel structure
#'
#' @description Creates a list of field names of a data frame with str
#' @family add-my-pet auxiliary functions
#' @param data data frame with fields
#' @param str string with field name
#' @return list of strings with fields including the field name in str
#' @examples nmpsd <- fieldnm_wtxt(data, "psd")
#' @export
fieldnm_wtxt <- function(data = list(), str = ""){

  nmaux <- names(data);
  fullnmaux <- list();
  fullnm <- list();

  while(length(nmaux) > 0) {
    if(nmaux[1] == str) {
      fullnm <- list(str);
    } else if(is.list(eval(parse(text = paste("data$", nmaux[1], sep = "")))))
      fullnmaux <- c(fullnmaux, nmaux[1]);
    if(length(nmaux) > 1)
      nmaux[1] <- NULL
    else
      nmaux <- list()
  }

  while(length(fullnmaux) > 0) {
    nmaux <- names(eval(parse(text = paste("data$", fullnmaux[1], sep = ""))));

    for(currentNm in nmaux) {
      if(currentNm == str) {
        fullnm <- c(fullnm, paste(fullnmaux[1], "$", str, sep = ""));
      } else if(is.list(eval(parse(text = paste("data$", fullnmaux[1], "$", currentNm, sep = "")))))
        fullnmaux <- c(fullnmaux, paste(fullnmaux[1], "$", str, sep = ""));
    }
    if(length(fullnmaux) > 1)
      fullnmaux[1] <- NULL
    else
      fullnmaux <- list()
  }

  return(fullnm)
}
