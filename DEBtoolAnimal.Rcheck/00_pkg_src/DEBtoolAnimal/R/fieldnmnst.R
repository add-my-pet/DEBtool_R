#' Creates a list of field names of a structure
#'
#' @description Creates a list of field names of a structure
#' @family add-my-pet auxiliary functions
#' @param st data frame with fields
#' @return list of list of strings with fields including the field name in str
#' @examples nst <- fieldnmnst(st)
#' @export
fieldnmnst <- function(st){

  nm <- as.list(names(st));
  baux <- 1;

  while(baux <= length(nm)) {
    for(currentField in nm)
      if(is.list(st[[currentField]])) {
        vaux <- names(st[[currentField]]);
        for(vauxField in vaux)
          nm[[length(nm)+1]] <- c(currentField, vauxField)
        #nm <- lapply(nm, function(x) x[sapply(x, length) >= 3]);
        nm <- nm[nm != currentField];
      } else
        baux <- baux + 1;
  }


  return(nm)
}
