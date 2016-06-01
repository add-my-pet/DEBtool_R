#' Sets automatically the weights for the data (to be used in a regression)
#'
#' @description computes weights for given data and adds it to the weight structure
#' @family add-my-pet auxiliary functions
#' @param data structure with data
#' @param weights structure with weights
#' @return structure with weights
#' @details computes weights for given data and adds it to the weight structure
#' for the zero-variate data y, the weight will be
#' \deqn{min(100, 1/ max(10^-6, y) ^2 \right)}{\min\left(100, \frac{1}{\max\left(10^{-6}, y\right ) ^2}}
#' for the uni-variate data y, the weight will be
#' \deqn{1/ N \bar{y}^2}{\frac{1}{N \bar{y}^2}}
#' @examples setweights(data)
#' @export
setweights <- function(data, weights = list()){
  nm <- names(data); # vector of cells with names of data sets

  for (dti in nm){
    dtSet = get(dti, data);
    if (is.null(weights[[dti]])) {
      if (length(dtSet) == 1) {
        weights[[dti]] <- min(100, 1/ max(1e-6, dtSet)^ 2);
      } else {
        ndt = dim(dtSet);
        meanval = mean(dtSet[,2]);
        weights[[dti]] <- rep(1/ meanval^2/ ndt[1], ndt[1]);
      }
    }
  }

  return(weights)

}
