#' Concatenates mydata files for several species
#'
#' @description Concatenates mydata files for several species
#' @family add-my-pet functions
#' @return structure with data, auxData, metaData, txtData and weights for several pets
#' @examples mydata_pets()
#' @export
mydata_pets <- function(){

  for(currentPet in pets) {
    list[data[[currentPet]], auxData[[currentPet]], metaData[[currentPet]], txtData[[currentPet]], weights[[currentPet]]] <- eval(parse(text = paste("mydata_ ", currentPet, "()", sep = "")))
  }

  return(list(data = data, auxData = auxData, metaData = metaData, txtData = txtData, weights = weights))

}
