#' Concatenates mydata files for several species
#'
#' @description Concatenates mydata files for several species
#' @family add-my-pet functions
#' @return structure with data, auxData, metaData, txtData and weights for several pets
#' @examples mydata_pets()
#' @export
mydata_pets <- function(){

  data <- list(); auxData <- list(); metaData <- list(); txtData <- list(); weights <- list();
  for(currentPet in pets) {
    list[data[[currentPet]], auxData[[currentPet]], metaData[[currentPet]], txtData[[currentPet]], weights[[currentPet]]] <- eval(parse(text = paste("mydata_", currentPet, "()", sep = "")))
  }



  return(list(data = data, auxData = auxData, metaData = metaData, txtData = txtData, weights = weights))
}

# if(pseudodata_pets == 1) {
#   ## remove pseudodata from species struture
#   data <- rmpseudodata(data);
#   txtData <- rmpseudodata(txtData);
#   weights <- rmpseudodata(weights);
#
#   ## set pseudodata and respective weights
#   list[dataTemp, unitsTemp, labelTemp, weightTemp] = addpseudodata();
#   data$psd <- dataTemp$psd;
#   weights$psd <- weightTemp$psd;
#   txtData$psd$units <- unitsTemp$psd;
#   txtData$psd$label <- labelTemp$psd;
# }
